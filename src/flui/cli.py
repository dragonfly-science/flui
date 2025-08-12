from datetime import datetime
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Annotated, TypeAlias, Union

import platformdirs
import typer
from loguru import logger
from rich import print
from rich.console import Console
from rich.progress import Progress
from typer import Option, Typer

from .app import FluiApp
from .dna import SegmentType, iter_reads, open_as_text
from .settings import Settings, SettingsError, get_settings
from .subtype import BarcodeSet, KmerIndexFlu, KmerSet

try:
    VERSION = version("flui-tui")
except PackageNotFoundError:
    # Fallback for development mode
    VERSION = "dev"


# https://github.com/Textualize/rich/issues/2416#issuecomment-1193773381
def _log_formatter(record) -> str:
    """Log message formatter."""
    color_map = {
        "TRACE": "dim blue",
        "DEBUG": "cyan",
        "INFO": "bold",
        "SUCCESS": "bold green",
        "WARNING": "yellow",
        "ERROR": "bold red",
        "CRITICAL": "bold white on red",
    }
    lvl_color = color_map.get(record["level"].name, "cyan")
    return (
        "[not bold green]{time:YYYY/MM/DD HH:mm:ss}[/not bold green]|{level.icon}"
        f"|[{lvl_color}]{{message}}[/{lvl_color}]"
    )


def init_logging(debug: bool):
    logger.remove()
    logdir = Path(platformdirs.user_log_dir("flui", ensure_exists=True))
    level = "DEBUG" if debug else "WARNING"
    logger.add(logdir / "flui.log", rotation="100 MB", level=level)

    console = Console()
    logger.add(
        console.print,
        level="WARNING",
        format=_log_formatter,  # type: ignore
    )


init_logging(False)


flui_app = Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
)

FastQDirOpt: TypeAlias = Annotated[
    Path, Option(exists=True, dir_okay=True, file_okay=False, help="FastQ folder")
]
FastaFileOpt: TypeAlias = Annotated[
    Path,
    Option(exists=True, dir_okay=False, file_okay=True, help="A fasta file"),
]
ExportCSVOpt: TypeAlias = Annotated[
    bool,
    Option(help="Export a CSV file"),
]
WorkerOpt: TypeAlias = Annotated[Union[int, None], Option(min=1, max=5)]


@flui_app.command()
def ui(
    run: FastQDirOpt,
    ref: FastaFileOpt,
    dump: ExportCSVOpt = True,
    workers: WorkerOpt = None,
):
    """Run the subtyping UI."""
    try:
        settings = get_settings()
    except SettingsError as e:
        print("Configuration file errors:")
        for error in e.readable_errors:
            print(f" - {error}")
        raise typer.Exit(code=1)  # noqa: B904

    if workers is not None:
        settings.workers = workers

    barcodes = BarcodeSet.create(root=run, ref_path=ref, ha_size=17, na_size=13)
    app = FluiApp(barcodes, settings)
    app.run()

    # Get the current date and time in a specific format
    now = datetime.now()
    str_now = now.strftime("%Y-%m-%d_%H-%M-%S")
    filename = f"flui_barcodes_{str_now}"
    pth = Path() / filename
    csv_path = pth.with_suffix(".csv")

    if dump:
        print(f"Writing CSV summary to [blue]`{csv_path}`[/blue]...")
        barcodes.write_csv_summary(pth.with_suffix(".csv"))

        json_pth = pth.with_suffix(".json")
        print(f"Write JSON overview to [blue]`{json_pth}[/blue]...")
        barcodes.write_json_summary(now, json_pth)

    print("Shutting down processes...")


# Prepend the version
ui.__doc__ = f"""
    [bold]FLUI v{VERSION}[/bold]

     {ui.__doc__}
"""

subtype_app = Typer(
    add_completion=False,
    no_args_is_help=True,
    rich_markup_mode="rich",
)

KmerSizeOpt: TypeAlias = Annotated[int, Option(min=7, max=31)]

LOCAL_SETTINGS = Settings()


@subtype_app.command()
def subtype(
    ref: FastaFileOpt,
    barcode: FastQDirOpt,
    ha_size: int = LOCAL_SETTINGS.ha_kmer_size,
    na_size: int = LOCAL_SETTINGS.na_kmer_size,
    max_files: int | None = None,
):
    """Find the subtype for a single barcode folder."""
    ref = ref.resolve()
    barcode = barcode.resolve()
    print(f"Using reference at `{ref}`")
    print(f"Reading FASTQ from folder `{barcode}`")

    ha_idx = KmerIndexFlu(location=ref, segment=SegmentType.HA, size=ha_size)
    ha_idx.load()
    na_idx = KmerIndexFlu(location=ref, segment=SegmentType.NA, size=na_size)
    na_idx.load()

    all_kmers = ha_idx.all_kmers | na_idx.all_kmers

    ha_reads = KmerSet(size=ha_size)
    na_reads = KmerSet(size=na_size)
    fastq_files = list(barcode.glob("*.fastq.gz"))

    if max_files is not None:
        fastq_files = fastq_files[:max_files]

    with Progress() as progress:
        task = progress.add_task("Reading FASTQ files...", total=len(fastq_files))
        for fastq_file in fastq_files:
            with open_as_text(fastq_file) as fd:
                for read in iter_reads(fd):
                    ha_reads.add_read(read.encode("ascii"), only=all_kmers)
                    na_reads.add_read(read.encode("ascii"), only=all_kmers)
            progress.update(task, advance=1)

    print(ha_idx.calc_match(ha_reads))
    print(na_idx.calc_match(na_reads))
