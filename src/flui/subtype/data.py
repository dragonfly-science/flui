import csv
import re
from datetime import datetime
from importlib.metadata import version
from pathlib import Path
from typing import ClassVar, Self

from loguru import logger
from pydantic import BaseModel, Field, computed_field

from flui.dna import SegmentType
from flui.settings import Settings, get_settings

from .kmer_index import KmerIndexFlu, SegmentMatch
from .kmer_set import KmerSet


class Reads(BaseModel):
    """A FastQ file containing a set of reads.

    For now, we assume that the file is named according to the RE below.
    """

    # The FASTQ file name is expected to be in this format:
    # `PAW31929_pass_barcode01_d026da3b_d8ed7bbf_0.fastq.gz`
    # TODO: it may be more reliable to read the first line from the file?

    re_fastq_name: ClassVar[re.Pattern] = re.compile(
        r"""
        ^(?P<flow_cell_id>[^_]+)_
        (?P<status>pass|fail)_
        barcode(?P<barcode_id>\d{2})_
        (?P<xxx_id>[a-f0-9]{8})?_?
        (?P<run_id>[a-f0-9]{8})_
        (?P<file_index>\d+)
        \.fastq\.gz$
        """,
        re.VERBOSE,
    )

    # These are taken from the files.
    # They should be the same for all files in a set (We don't check).
    flow_cell_id: str
    barcode_id: str
    run_id: str

    # Extracted from the filename.
    file_index: int

    # The path to the file
    path: Path = Field(repr=False)

    # The path relative to the container the file
    rel_path: Path | None = Field(repr=False, default=None)

    ha_kmers: KmerSet | None = None
    na_kmers: KmerSet | None = None

    count: int = 0

    @property
    def key(self) -> str:
        return "/".join([self.run_id, self.barcode_id, str(self.file_index)])

    @property
    def barcode_key(self) -> str:
        return f"{self.run_id}/{self.barcode_id}"

    @classmethod
    def from_path(cls, path: Path) -> Self | None:
        """Create a Reads from a fastq file."""
        match = cls.re_fastq_name.match(path.name)
        if not match:
            logger.warning(f"Filename does not match expectations: `{path.name}`")
            return None

        slf = cls(
            flow_cell_id=match.group("flow_cell_id"),
            barcode_id=match.group("barcode_id"),
            file_index=int(match.group("file_index")),
            run_id=f"{match.group('run_id')}",
            path=path,
        )
        logger.info(f"loaded {slf}")
        return slf


class Barcode(BaseModel):
    """A barcode folder, containing a set of reads."""

    folder: Path
    run_id: str
    barcode_id: str
    flow_cell_id: str

    # Don't show this.
    reads: dict[int, Reads] = Field(repr=False, default_factory=dict)
    ha_kmers: KmerSet = Field(repr=False)
    na_kmers: KmerSet = Field(repr=False)

    # Summarized information that we update over time.
    reads_count: int = 0
    ha_kmer_count: int = 0
    na_kmer_count: int = 0

    # Our Matches. These will be updated as we add reads.
    ha_match: SegmentMatch | None = None
    na_match: SegmentMatch | None = None

    def get_assigned(self, segment_type: SegmentType, settings: Settings) -> str:
        if segment_type is SegmentType.HA:
            sm = self.ha_match
            sc = self.ha_kmer_count
        else:
            sm = self.na_match
            sc = self.na_kmer_count

        if sc < settings.minimum_kmers or sm is None:
            return "?"
        return sm.assigned_subtype(settings.minimum_score, settings.minimum_gap)

    @computed_field
    @property
    def ha_subtype(self) -> str:
        return self.ha_match.best_subtype() if self.ha_match else "?"

    @computed_field
    @property
    def na_subtype(self) -> str:
        return self.na_match.best_subtype() if self.na_match else "?"

    @computed_field
    @property
    def ha_score(self) -> float:
        return self.ha_match.best_score() if self.ha_match else 0.0

    @computed_field
    @property
    def na_score(self) -> float:
        return self.na_match.best_score() if self.na_match else 0.0

    @computed_field
    @property
    def key(self) -> str:
        return f"{self.run_id}/{self.barcode_id}"

    @computed_field
    @property
    def file_count(self) -> int:
        return len(self.reads)

    @computed_field
    @property
    def done_count(self) -> int:
        return sum(1 for r in self.reads.values() if r.ha_kmers is not None)

    def add_reads(self, reads: Reads) -> bool:
        if reads.file_index in self.reads:
            logger.warning(f"Duplicate reads index `{reads.file_index}`")
            return False

        self.reads[reads.file_index] = reads
        return True

    def update_kmers_from_reads(self, reads: Reads):
        if reads.file_index not in self.reads:
            logger.error(f"Cannot find read index `{reads.file_index}`")
            return

        if reads.ha_kmers is None or reads.na_kmers is None:
            logger.error("Kmers not calculated!")
            return

        # We just need to recalculate the barcode info.
        self.ha_kmers.add_set(reads.ha_kmers)
        self.na_kmers.add_set(reads.na_kmers)

        # Keep some running totals as we read in.
        self.reads_count += reads.count
        self.ha_kmer_count += reads.ha_kmers.distinct
        self.na_kmer_count += reads.na_kmers.distinct


class BarcodeSet(BaseModel):
    """A set of barcode runs from a sequencer.

    This will scan everything within the root folder.
    """

    # Where it is based.
    root: Path
    ha_index: KmerIndexFlu
    na_index: KmerIndexFlu

    barcodes: dict[str, Barcode] = Field(repr=False, default_factory=dict)

    @classmethod
    def create(cls, root: Path, ref_path: Path, ha_size: int, na_size: int) -> Self:
        ha_index = KmerIndexFlu(location=ref_path, segment=SegmentType.HA, size=ha_size)
        na_index = KmerIndexFlu(location=ref_path, segment=SegmentType.NA, size=na_size)
        return cls(root=root, ha_index=ha_index, na_index=na_index)

    def add_reads(self, reads: Reads) -> Barcode | None:
        # Sanity checking - resolve symlinks for proper path comparison
        try:
            resolved_root = self.root.resolve()
            resolved_reads_path = reads.path.resolve()
            rel_path = resolved_reads_path.relative_to(resolved_root)
        except ValueError:
            logger.warning(
                f"Reads path not in run folder: `{reads.path}` (root: {self.root})"
            )
            return None

        # Update this.
        reads.rel_path = rel_path

        # Create a barcode if we don't have one.
        barcode = self.barcodes.get(reads.barcode_key)

        if barcode is None:
            barcode = Barcode(
                barcode_id=reads.barcode_id,
                folder=reads.path.parent,
                flow_cell_id=reads.flow_cell_id,
                run_id=reads.run_id,
                ha_kmers=KmerSet(size=self.ha_index.size),
                na_kmers=KmerSet(size=self.na_index.size),
            )
            self.barcodes[barcode.key] = barcode

        # Now add the reads.
        if barcode.add_reads(reads):
            return barcode

        return None

    def update_reads(self, reads: Reads) -> Barcode:
        """Force everything to recalculate, related to these reads."""
        # Locate the reaod
        barcode = self.barcodes[reads.barcode_key]
        barcode.update_kmers_from_reads(reads)
        return barcode

        # recalculate the subtype and confidence.
        # barcode.subtype = self.index.get_subtype(barcode.kmers)

    def write_csv_summary(self, path: Path):
        st = get_settings()
        sorted_barcodes = sorted(self.barcodes.values(), key=lambda b: b.key)
        with path.open("w", encoding="utf-8") as fd:
            writer = csv.writer(fd)
            writer.writerow(
                [
                    "run_id",
                    "flow_cell_id",
                    "barcode_id",
                    "file_count",
                    "done_count",
                    "ha_kmer_count",
                    "ha_subtype",
                    "ha_score",
                    "na_kmer_count",
                    "na_subtype",
                    "na_score",
                ]
            )

            for bc in sorted_barcodes:
                writer.writerow(
                    [
                        bc.run_id,
                        bc.flow_cell_id,
                        bc.barcode_id,
                        bc.file_count,
                        bc.done_count,
                        bc.ha_kmer_count,
                        bc.get_assigned(SegmentType.HA, st),
                        bc.ha_score,
                        bc.na_kmer_count,
                        bc.get_assigned(SegmentType.NA, st),
                        bc.na_score,
                    ]
                )

    def write_json_summary(self, when: datetime, path: Path):
        """Capture an overview of the barcode set."""
        j = JsonSummary.from_barcode_set(when, self)
        with path.open("w", encoding="utf-8") as f:
            f.write(j.model_dump_json(indent=2))


class JsonSummary(BaseModel):
    """A summary of the barcode set."""

    flui_version: str
    when: datetime
    where: str

    @classmethod
    def from_barcode_set(cls, when: datetime, barcode: BarcodeSet) -> Self:
        where = barcode.root.as_posix()
        return cls(where=where, when=when, flui_version=version("flui-tui"))
