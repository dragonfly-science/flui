import argparse
import asyncio
import os
import random
import shutil
from pathlib import Path


async def simseq(
    src: Path,
    dst: Path,
    sleep: float = 1.0,
):
    """Simulate incoming reads from a reference."""
    src = src.resolve()
    dst = dst.resolve()

    all_reads = list(src.glob("**/*.fastq.gz"))
    if not all_reads:
        print(f"No FASTQ files found in {src}")
        return

    print(f"Found {len(all_reads)} FASTQ files to simulate")
    random.shuffle(all_reads)

    for fastq in all_reads:
        wait = random.expovariate(1 / sleep)
        print(f"Waiting {wait:.2f}s before copying {fastq.name}")
        await asyncio.sleep(wait)

        relative_path = fastq.relative_to(src)
        dest = dst / relative_path
        dest.parent.mkdir(parents=True, exist_ok=True)
        # Copy the file to the destination
        shutil.copy2(fastq, dest)
        print(f"Copied: {fastq} -> {dest}")


def clean_destination(dst: Path) -> bool:
    """Clean the destination directory after user confirmation."""
    if not dst.exists():
        return True

    files = list(dst.rglob("*"))
    if not files:
        return True

    print(f"Destination directory {dst} contains {len(files)} files/directories.")
    response = input("Remove all contents? [y/N]: ").strip().lower()

    if response in {"y", "yes"}:
        print(f"Cleaning {dst}...")
        for item in dst.iterdir():
            if item.is_dir():
                shutil.rmtree(item)
            else:
                item.unlink()
        print("Destination cleaned.")
        return True
    else:
        print("Clean cancelled.")
        return False


def load_env_file(env_file: Path = Path(".env")) -> dict:
    """Load environment variables from a .env file."""
    env_vars = {}
    if env_file.exists():
        with env_file.open(encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    key, _, value = line.partition("=")
                    env_vars[key.strip()] = value.strip()
    return env_vars


def main():
    # Load .env file
    env_vars = load_env_file()

    # Get defaults from environment or .env file
    default_src = env_vars.get("SIMSEQ_SRC") or os.getenv("SIMSEQ_SRC")
    default_dst = env_vars.get("SIMSEQ_DST") or os.getenv("SIMSEQ_DST")
    default_sleep = float(
        env_vars.get("SIMSEQ_SLEEP") or os.getenv("SIMSEQ_SLEEP", "1.0")
    )

    parser = argparse.ArgumentParser(
        description=(
            "Simulate incoming sequencing files by copying FASTQ files "
            "with random delays"
        )
    )
    parser.add_argument(
        "src",
        type=Path,
        nargs="?" if default_src else None,
        default=Path(default_src) if default_src else None,
        help=(
            f"Source directory containing FASTQ files"
            f"{' (default from .env)' if default_src else ''}"
        ),
    )
    parser.add_argument(
        "dst",
        type=Path,
        nargs="?" if default_dst else None,
        default=Path(default_dst) if default_dst else None,
        help=(
            f"Destination directory to copy files to"
            f"{' (default from .env)' if default_dst else ''}"
        ),
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=default_sleep,
        help=(
            f"Average sleep time between file copies "
            f"(seconds, default: {default_sleep})"
        ),
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Remove all contents from destination directory before starting",
    )

    args = parser.parse_args()

    if not args.src:
        parser.error("src is required (provide as argument or set SIMSEQ_SRC in .env)")
    if not args.dst:
        parser.error("dst is required (provide as argument or set SIMSEQ_DST in .env)")

    if not args.src.exists():
        print(f"Error: Source directory {args.src} does not exist")
        return 1

    if not args.src.is_dir():
        print(f"Error: Source {args.src} is not a directory")
        return 1

    # Clean destination if requested
    if args.clean and not clean_destination(args.dst):
        print("Exiting due to cancelled clean operation.")
        return 1

    print(f"Simulating sequencing from {args.src} to {args.dst}")
    print(f"Average delay: {args.sleep}s")

    try:
        asyncio.run(simseq(args.src, args.dst, args.sleep))
    except KeyboardInterrupt:
        print("\nSimulation interrupted")
        return 1

    print("Simulation complete")
    return 0


if __name__ == "__main__":
    exit(main())
