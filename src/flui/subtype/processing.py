import asyncio
import random
from asyncio import Queue, Task
from concurrent.futures import Future
from concurrent.futures import ProcessPoolExecutor as Executor
from datetime import datetime
from enum import Enum, auto
from pathlib import Path

from loguru import logger
from pydantic import BaseModel, Field
from watchfiles import Change, awatch

from flui.dna import iter_reads, open_as_text

from .data import Barcode, BarcodeSet, Reads
from .kmer_index import KmerIndexFlu, SegmentMatch
from .kmer_set import KmerSet

# This is our remote data store.
WORKER_DATA: list[KmerIndexFlu] = []
HA_INDEX = 0
NA_INDEX = 1


def _initialize_kmers(ha_index: KmerIndexFlu, na_index: KmerIndexFlu):
    # We load the indexes in the worker processes.
    ha_index.load()
    na_index.load()
    WORKER_DATA.append(ha_index)
    WORKER_DATA.append(na_index)


# Break this out to make it play nicely with the ProcessPoolExecutor
def _process_reads(fastq: Reads) -> tuple[int, KmerSet, KmerSet]:
    if not WORKER_DATA:
        msg = "KmerIndexes not initialised"
        raise RuntimeError(msg)
    ha_index = WORKER_DATA[HA_INDEX]
    na_index = WORKER_DATA[NA_INDEX]

    ha_reads = KmerSet(size=ha_index.size)
    na_reads = KmerSet(size=na_index.size)
    read_count = 0
    with open_as_text(fastq.path) as fd:
        for read in iter_reads(fd):
            read_count += 1
            ha_reads.add_read(read.encode("ascii"), only=ha_index.all_kmers)
            na_reads.add_read(read.encode("ascii"), only=na_index.all_kmers)
    return read_count, ha_reads, na_reads


def _match_barcode(
    ha_kmers: KmerSet, na_kmers: KmerSet
) -> tuple[SegmentMatch, SegmentMatch]:
    if not WORKER_DATA:
        msg = "KmerIndexes not initialised"
        raise RuntimeError(msg)
    ha_index = WORKER_DATA[HA_INDEX]
    na_index = WORKER_DATA[NA_INDEX]

    ha = ha_index.calc_match(ha_kmers)
    na = na_index.calc_match(na_kmers)
    return ha, na


def filter_fastq(change: Change, pth_str: str):
    pth = Path(pth_str)
    result = change == Change.added and pth.is_file() and pth_str.endswith(".fastq.gz")
    logger.debug(f"Filter check: {change} {pth_str} -> {result}")
    return result


class BarcodeUpdateKind(Enum):
    NEW_READS = auto()
    COMPLETE_READS = auto()
    UPDATE_MATCHES = auto()

    def as_text(self) -> str:
        return self.name.title().replace("_", " ")


class BarcodeUpdate(BaseModel):
    when: datetime = Field(default_factory=datetime.now)
    kind: BarcodeUpdateKind
    barcode: Barcode
    reads: Reads | None = None

    def __rich__(self) -> str:
        text = f"[dim][{self.when:%Y-%m-%d %H:%M:%S}][/dim] "
        text += f"[bold]{self.kind.as_text():<14}[/bold]: "
        text += self.barcode.key
        idx = self.reads.file_index if self.reads else -1

        match self.kind:
            case BarcodeUpdateKind.NEW_READS | BarcodeUpdateKind.COMPLETE_READS:
                text += f", file num {idx:02d}"

            case BarcodeUpdateKind.UPDATE_MATCHES:
                text += f", total reads {self.barcode.reads_count:,}"

        return text


class BarcodeProcessor:
    """Process a set of FastQ files."""

    def __init__(self, *, barcode_set: BarcodeSet, max_workers: int = 8):
        self.barcode_set = barcode_set
        self.read_executor = Executor(
            max_workers=max_workers,
            initializer=_initialize_kmers,
            initargs=(
                barcode_set.ha_index,
                barcode_set.na_index,
            ),
        )
        self.score_executor = Executor(
            max_workers=1,
            initializer=_initialize_kmers,
            initargs=(
                barcode_set.ha_index,
                barcode_set.na_index,
            ),
        )
        self.future_reads: dict[Future, Reads] = {}
        self.future_scoring: dict[Future, Barcode] = {}

        self.shutting_down = False

        # We'll boot these up.
        self.watch_task: Task | None = None
        self.scan_task: Task | None = None
        self.updates_q: Queue[BarcodeUpdate] = Queue()

    def start(self):
        # Start the watch before the scan, so we don't miss anything.
        # We should automatically detect if we're doubling up on the same file.
        logger.info(f"Starting file processor for {self.barcode_set.root}")
        self.watch_task = asyncio.create_task(self._watch())
        self.scan_task = asyncio.create_task(self._scan())

    async def get_updates(self) -> list[BarcodeUpdate]:
        updates = []
        while True:
            try:
                ret = self.updates_q.get_nowait()
                logger.info(f"Retrieved update from queue: {ret}")
            except asyncio.QueueEmpty:
                break
            updates.append(ret)
            await asyncio.sleep(0)
        if updates:
            logger.info(f"Returning {len(updates)} updates to TUI")
        return updates

    async def _scan(self):
        """Scan the barcode directory for existing FastQ files."""
        fastq_paths = list(self.barcode_set.root.glob("**/*.fastq.gz"))
        random.shuffle(fastq_paths)
        for pth in fastq_paths:
            if self.shutting_down:
                break
            await self._submit_new_reads(pth)

    async def _watch(self):
        """Use the watchfiles library to watch for new FastQ files."""
        logger.info(f"Starting watch on directory: {self.barcode_set.root}")

        # Check if directory exists before starting watch
        if not self.barcode_set.root.exists():
            logger.error(f"Watch directory does not exist: {self.barcode_set.root}")
            return

        if not self.barcode_set.root.is_dir():
            logger.error(f"Watch path is not a directory: {self.barcode_set.root}")
            return

        try:
            async for changes in awatch(
                self.barcode_set.root, watch_filter=filter_fastq
            ):
                if self.shutting_down:
                    logger.info("Watch shutting down")
                    break
                logger.info(f"Watch detected {len(changes)} changes: {changes}")
                for _, pth_str in changes:
                    logger.info(f"Processing new file: {pth_str}")
                    await self._submit_new_reads(Path(pth_str))
        except Exception as e:
            logger.error(f"Watch failed with error: {e}")
            raise

    async def _submit_new_reads(self, pth: Path):
        if self.shutting_down:
            logger.debug(f"Skipping {pth} - shutting down")
            return

        logger.debug(f"Attempting to create reads from {pth}")
        reads = Reads.from_path(pth)
        if reads:
            logger.info(f"Created reads from {pth}: count={reads.count}")
            barcode = self.barcode_set.add_reads(reads)
            if barcode:
                logger.info(f"Added reads to barcode {barcode.key}")
                update = BarcodeUpdate(
                    kind=BarcodeUpdateKind.NEW_READS, barcode=barcode, reads=reads
                )
                logger.info(f"Putting update in queue: {update}")
                await self.updates_q.put(update)
                # Submit them for processing.
                future = self.read_executor.submit(_process_reads, reads)
                self.future_reads[future] = reads
                future.add_done_callback(self._reads_done)
            else:
                logger.warning(f"Failed to add reads to barcode set for {pth}")
        else:
            logger.warning(f"Failed to create reads from {pth}")
        #     await self.updates_q.put((BarcodeUpdate.NEW_READS, None))

    def _reads_done(self, future: Future):
        if self.shutting_down or future.cancelled():
            return

        cnt, ha, na = future.result()
        reads = self.future_reads.pop(future)
        reads.count = cnt
        reads.ha_kmers = ha
        reads.na_kmers = na
        barcode = self.barcode_set.update_reads(reads)
        self.updates_q.put_nowait(
            BarcodeUpdate(
                kind=BarcodeUpdateKind.COMPLETE_READS, barcode=barcode, reads=reads
            )
        )
        # Rescore the barcode
        self._submit_matching(barcode)

    def _submit_matching(self, barcode: Barcode):
        if self.shutting_down:
            return

        # Already there? Well then, don't bother.
        if barcode in self.future_scoring.values():
            return

        future = self.score_executor.submit(
            _match_barcode, barcode.ha_kmers, barcode.na_kmers
        )
        self.future_scoring[future] = barcode
        future.add_done_callback(self._matching_done)

    def _matching_done(self, future: Future):
        if self.shutting_down or future.cancelled():
            return

        barcode = self.future_scoring.pop(future)
        ha, na = future.result()
        barcode.ha_match = ha
        barcode.na_match = na
        self.updates_q.put_nowait(
            BarcodeUpdate(kind=BarcodeUpdateKind.UPDATE_MATCHES, barcode=barcode)
        )

    def shutdown(self):
        """Shuts down the tasks and executors."""
        self.shutting_down = True
        for t in (self.watch_task, self.scan_task):
            if t is not None and not t.done():
                t.cancel()

        self.read_executor.shutdown(wait=False, cancel_futures=True)
        self.score_executor.shutdown(wait=False, cancel_futures=True)
