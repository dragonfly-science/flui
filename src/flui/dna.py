from __future__ import annotations

import gzip
import lzma
import string
from dataclasses import dataclass
from enum import IntEnum
from typing import TYPE_CHECKING, Self, TextIO, TypeAlias, cast

if TYPE_CHECKING:
    from collections.abc import Iterator
    from pathlib import Path

DNA_LETTERS = "ACGTWSMKRYBDHVN-"
DNA_LETTERS_BYTES = DNA_LETTERS.encode("ascii")
WHITE_TRANSLATION = str.maketrans("", "", string.whitespace)
ILLEGAL_TRANSLATION = str.maketrans("", "", DNA_LETTERS)


def get_illegal_from_normalized(s: str) -> str:
    """If you remove all the legal letters, you get this."""
    return s.translate(ILLEGAL_TRANSLATION)


def normalize_bases(s: str) -> str:
    normed = s.translate(WHITE_TRANSLATION).upper()

    # Error check here.
    bad_chars = get_illegal_from_normalized(normed)
    if bad_chars:
        msg = f"Invalid DNA sequence: it contains '{bad_chars}'"
        raise ValueError(msg)

    return normed


@dataclass(frozen=True, slots=True, repr=False)
class DNAString:
    """A simple DNA string, that has been normalized."""

    bases: str

    def __repr__(self) -> str:
        return f"DNAString({self.bases[:10]}...)"

    def to_bytes(self) -> bytes:
        return self.bases.encode("ascii")


FASTA_START = ">"
NameAndDNA: TypeAlias = tuple[str, DNAString]


class FastaSyntaxError(Exception):
    """Fasta file problems."""


def iter_fasta(fd: TextIO) -> Iterator[NameAndDNA]:
    # We prime everything by getting a description line.
    # Find the first non-blank line.
    line_number = 1
    while raw_line := fd.readline():
        line_number += 1
        desc_line = raw_line.strip()
        if desc_line:
            break
    else:
        msg = "Fasta file contains no sequences."
        raise FastaSyntaxError(msg)

    if desc_line[0] != FASTA_START:
        msg = "The first sequence does not begin with '>'"
        raise FastaSyntaxError(msg)

    # Loop over sequences
    while True:
        seq_line = ""
        lines = []

        # Loop over lines
        while raw_line := fd.readline():
            line_number += 1
            seq_line = raw_line.strip()
            if not seq_line:
                # Skip blank line. This goes beyond the standard.
                continue

            if seq_line[0] == FASTA_START:
                # We're starting a new sequence.
                break

            # Otherwise, append to the bases.
            try:
                normed = normalize_bases(seq_line)
            except ValueError as e:
                msg = (
                    f"Invalid chars in fasta file: {e}"
                    f" at line {line_number},"
                    f" while reading sequence `{desc_line[1:]}`."
                )
                raise FastaSyntaxError(msg) from e

            line_number += 1
            lines.append(normed)

        # We either ran out of data, or we found a new description line.
        # NOTE: we don't need to validate the DNAString -- we've done it above.
        yield desc_line[1:], DNAString("".join(lines))

        if not raw_line:
            # EOF
            break

        # Set up the description for the next sequence.
        desc_line = seq_line


def iter_reads(fd: TextIO) -> Iterator[str]:
    """Reading of FASTQ files full of reads."""
    while header := fd.readline():
        if header[0] != "@":
            msg = f"Expected a header line starting with '@', got {header}"
            raise ValueError(msg)

        # The sequence
        read = fd.readline().strip()

        # + line (typically empty).
        plus = fd.readline()
        if plus[0] != "+":
            msg = f"Expected a '+' line, got {plus}"
            raise ValueError(msg)

        # Quality line. We ignore for now.
        fd.readline()

        # Now we have checked have everything, yield the read.
        yield read


def open_as_text(path: Path, mode: str = "rt") -> TextIO:
    if "t" not in mode:
        # This makes the cast okay..
        msg = "Mode must be a text mode (include 't')"
        raise ValueError(msg)
    match path.suffix:
        case ".gz":
            return cast("TextIO", gzip.open(path, mode))
        case ".xz":
            return cast("TextIO", lzma.open(path, mode))
        case _:
            return cast("TextIO", path.open(mode))


class SegmentTypeError(Exception):
    pass


class SegmentType(IntEnum):
    """These names and numbers are related (not like a normal enum)."""

    PB2 = 1
    PB1 = 2
    PA = 3
    HA = 4
    NP = 5
    NA = 6
    MP = 7
    NS = 8

    def to_string(self):
        """Clearer, as I can never remember this."""
        return self.name

    @classmethod
    def parse(cls, text: str) -> Self:
        """Support parsing of integers or strings."""
        text = text.strip().upper()

        # Try a numeric value first.
        try:
            v = int(text)
        except ValueError:
            pass
        else:
            # It should be a valid value
            try:
                return cls(v)
            except ValueError as e:
                msg = f"{v} is not a valid segment"
                raise SegmentTypeError(msg) from e

        try:
            return cls[text]
        except KeyError as e:
            msg = f"{text} is not a valid segment"
            raise SegmentTypeError(msg) from e
