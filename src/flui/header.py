import re
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, ClassVar, Self

from loguru import logger
from pydantic import BaseModel, field_validator

from .dna import SegmentType, SegmentTypeError, iter_fasta, open_as_text

RE_CLASSIFY_SUBTYPE = re.compile(r"H(?P<ha>\d+)N(?P<na>\d+)")

RE_SEGMENT_TYPES = "|".join(s.to_string() for s in SegmentType)

# Various regexes for parsing parts of the header.
# These will typically be separated with |.
# We use these to build different header recognition regexes.

# A simple UUID recogniser (we don't check lengths).
RE_UUID = r"\s*(?P<uuid>[a-z0-9-]+)\s*"

# Subtypes can be empty. Otherwise, we except numbers and uppercase.
# Inspection reveals a few other chars like ',/?'.
RE_SUBTYPE = r"(?P<subtype>[A-Z0-9?,/_]{0,9}\s*)"

# Segment number is a single integer. We strip whitespace.
RE_SEGMENT_NUM = r"\s*(?P<segment_num>\d{0,1})\s*"

# Segment name is one of the types from the Enum.
RE_SEGMENT = rf"\s*(?P<segment>{RE_SEGMENT_TYPES})\s*"

# Mostly, virus names should have some forward /s in them somewhere.
RE_SLASH_IN_NAME = r"\s*(?P<name>[^|]+?/[^|]+)\s*"

# But NCBI sometimes breaks this rule.
# Typically, it has a string at the beginning, in whatever case.
RE_NCBI_NAME = r"\s*(?i:Influenza A virus)\s*(?P<name>[^|]+?)\s*"

# The NZ names have a bunch of information we can extract.
RE_TYPE = r"(?P<type>[A-Za-z])"
RE_HOST = r"(?P<host>[A-Za-z _-]+)"
RE_COUNTRY = r"(?P<country>[A-Za-z _-]+)"  # NZL or NZ, for now.
RE_ISOLATE = r"(?P<isolate>[A-Za-z0-9_-]+)"
RE_YEAR = r"(?P<year>\d{4})"
RE_NZ_NAME = rf"{RE_TYPE}/{RE_HOST}/{RE_COUNTRY}/{RE_ISOLATE}/{RE_YEAR}\({RE_SUBTYPE}\)"

# This should be an NCBI or GISAID key.
# These should be in the format EPINNNNN or a ABNNNN.N
# We set maximum size of 15 (currently 10).
RE_ACCESSION = r"\s*(?P<accession>[A-Z0-9._]{0,20})\s*"
RE_GISAID_ISOLATE = r"\s*(?P<isolate>[A-Z0-9_]{0,20})\s*"


class FastaHeaderError(Exception):
    """Fasta header is in unexpected format."""


class FastaHeader(BaseModel, ABC):
    # A unique identifier for the record.
    sequence_key: str

    # Some name (maybe not unique as it hold across segments)
    name: str

    # A SegmentType, if present.
    # We can parse this from an int or the string (e.g. "HA")
    segment: SegmentType | None

    # Free form subtype. We use None (rather than "") to make nulls in the database.
    subtype: str | None

    # These are defaulted to None, as they are not always present.
    isolate_key: str | None = None
    location: str | None = None
    collection_date: str | None = None
    host: str | None = None

    # Registry of all Fasta header types.
    registry: ClassVar[dict[str, type["FastaHeader"]]] = {}

    @classmethod
    def samples(cls) -> list[str]:
        return []

    @staticmethod
    def detect(header: str) -> str:
        for nm, cls in FastaHeader.registry.items():
            try:
                cls.parse(header)
                return nm
            except FastaHeaderError:
                continue

        msg = f"No match for: {header}"
        raise FastaHeaderError(msg)

    @staticmethod
    def get_parser(header: str) -> type["FastaHeader"] | None:
        try:
            nm = FastaHeader.detect(header)
        except FastaHeaderError:
            return None

        return FastaHeader.registry[nm]

    @staticmethod
    def get_parser_for_file(path: Path) -> type["FastaHeader"] | None:
        with open_as_text(path) as fd:
            for nm, _ in iter_fasta(fd):
                return FastaHeader.get_parser(nm)
        return None

    @classmethod
    def get_name(cls) -> str:
        return cls.__name__.replace("Header", "")

    @classmethod
    def __pydantic_init_subclass__(cls, **kwargs: Any) -> None:
        """This fills out the registry."""
        super().__init_subclass__(**kwargs)
        cls.registry[cls.get_name()] = cls

    @classmethod
    def do_parse(cls, header_re: re.Pattern, header: str) -> re.Match:
        match = header_re.match(header)

        if not match:
            msg = f"Expected {cls.__name__} type, but got non-matching: `{header}`"
            raise FastaHeaderError(msg)
        return match

    @classmethod
    @abstractmethod
    def parse(cls, header: str) -> Self: ...

    @field_validator("host", mode="before")
    @classmethod
    def to_lowercase(cls, raw: str | None) -> str | None:
        if raw is None:
            return None
        return raw.lower()

    @field_validator("subtype", mode="before")
    @classmethod
    def clean_segment(cls, v: str) -> str | None:
        if v == "":
            return None
        return v.upper()

    def has_valid_subtype(self) -> bool:
        return (
            self.subtype is not None
            and RE_CLASSIFY_SUBTYPE.fullmatch(self.subtype) is not None
        )

    def classify_segment(self) -> str | None:
        """Classify the header into segment/subtype."""
        if self.subtype is None:
            logger.warning("Attempt to classify segment with blank subtype")
            return None

        subtype_match = RE_CLASSIFY_SUBTYPE.fullmatch(self.subtype)
        if subtype_match is None:
            logger.warning(f"`{self.subtype}` does not match specify H and N types")
            return None

        # We group the sequences by segment and appropriate subtype.
        ha = subtype_match.group("ha")
        na = subtype_match.group("na")

        if self.segment == SegmentType.HA:
            segment_type = f"H{ha}"
        elif self.segment == SegmentType.NA:
            segment_type = f"N{na}"
        else:
            logger.warning(
                "Attempt to classify segment which is neither HA or NA:"
                f"(it is segment `{self.segment}`)"
            )
            segment_type = None

        return segment_type


class NZHeader(FastaHeader):
    """NZ Fasta Headers."""

    @classmethod
    def samples(cls) -> list[str]:
        return [
            "1|PB2|A/Mallard/NZL/W19_10_79/2019(H3N8)",
            "1|PB2|A/Mallard/NZ/W11_214_c152/2011(H5N2)",
        ]

    HEADER_RE: ClassVar[re.Pattern] = re.compile(
        r"\|".join(  # noqa: FLY002
            [
                RE_SEGMENT_NUM,
                RE_SEGMENT,
                RE_NZ_NAME,
            ]
        )
    )

    @classmethod
    def parse(cls, header: str) -> Self:
        match = cls.do_parse(cls.HEADER_RE, header)

        subtype = match.group("subtype")
        parts = {
            nm: match.group(nm) for nm in ["type", "host", "country", "isolate", "year"]
        }
        name = "/".join(parts.values()) + f"({subtype})"

        segment = SegmentType.parse(match.group("segment"))
        key = f"{parts['isolate']}/{segment.value}"
        return cls(
            segment=segment,
            name=name,
            host=parts["host"],
            location=parts["country"],
            sequence_key=key,
            isolate_key=parts["isolate"],
            subtype=subtype,
            collection_date=parts["year"],
        )


class NCBIHeader(FastaHeader):
    """NCBI Fasta Headers."""

    @classmethod
    def samples(cls) -> list[str]:
        return [
            (
                " |Influenza A virus (A/goose/Guangdong/1/1996(H5N1)) polymerase (PA) "
                "and PA-X protein (PA-X) genes, complete cds|NC_007359.1|H5N1"
            ),
            (
                "2 |Influenza A virus (A/goose/Guangdong/1/1996(H5N1)) polymerase "
                "(PB1) and PB1-F2 protein (PB1-F2) genes, complete cds|NC_007358.1|H5N1"
            ),
        ]

    HEADER_RE: ClassVar[re.Pattern] = re.compile(
        r"\|".join(  # noqa: FLY002
            [
                RE_SEGMENT_NUM,
                RE_NCBI_NAME,
                RE_ACCESSION,
                RE_SUBTYPE,
            ]
        )
    )

    @classmethod
    def parse(cls, header: str) -> Self:
        match = cls.do_parse(cls.HEADER_RE, header)
        try:
            segment = SegmentType.parse(match.group("segment_num"))
        except SegmentTypeError:
            segment = None

        return cls(
            segment=segment,
            sequence_key=match.group("accession"),
            name=match.group("name"),
            subtype=match.group("subtype"),
        )


RE_TRAILING_NONDIGITS = re.compile(r"\D+$")


class GISAIDHeader(FastaHeader):
    """GISAID Headers."""

    @classmethod
    def samples(cls) -> list[str]:
        return [
            "2|PB1|EPI_ISL_131202|A/duck/Alberta/35/1976|CY009609|A/H3N",
            "5|NP|EPI_ISL_19136018|A/gull/France/23P003123/2023|EPI3283011|A_/_H5N1",
        ]

    HEADER_RE: ClassVar[re.Pattern] = re.compile(
        r"\|".join(  # noqa: FLY002
            [
                RE_SEGMENT_NUM,
                RE_SEGMENT,
                RE_GISAID_ISOLATE,
                RE_SLASH_IN_NAME,
                RE_ACCESSION,
                RE_SUBTYPE,
            ]
        )
    )

    @classmethod
    def parse(cls, header: str) -> Self:
        match = cls.do_parse(cls.HEADER_RE, header)

        # We need to remove type, as gisaid uses (A/H3N2) for the "type"
        subtype = match.group("subtype")
        if "/" in subtype:
            subtype = subtype.split("/")[-1]
            subtype = subtype.replace("_", "")
        else:
            subtype = ""

        isolate = match.group("isolate")
        segment = SegmentType.parse(match.group("segment"))

        # deal with the name
        name = match.group("name")
        parts = name.split("/")

        host = parts[1] if len(parts) > 2 else None

        # Let's try and pull out a date
        # Grab last component. Remove any brackets (sometimes contains subtype).
        date = parts[-1].strip()
        date = RE_TRAILING_NONDIGITS.sub("", date)
        try:
            int(date)
        except ValueError:
            date = None

        key = match.group("accession")
        return cls(
            segment=segment,
            name=name,
            host=host,
            sequence_key=key,
            subtype=subtype,
            isolate_key=isolate,
            collection_date=date,
        )


class SimpleGISAIDHeader(FastaHeader):
    """Original GISAID Header."""

    @classmethod
    def samples(cls) -> list[str]:
        return [
            "1|PB2|EPI_ISL_131202|H4N6",
        ]

    HEADER_RE: ClassVar[re.Pattern] = re.compile(
        r"\|".join(  # noqa: FLY002
            [
                RE_SEGMENT_NUM,
                RE_SEGMENT,
                RE_GISAID_ISOLATE,
                RE_SUBTYPE,
            ]
        )
    )

    @classmethod
    def parse(cls, header: str) -> Self:
        match = cls.do_parse(cls.HEADER_RE, header)
        isolate = match.group("isolate")
        segment = SegmentType.parse(match.group("segment"))

        key = f"{isolate}/{segment.value}"
        return cls(
            segment=segment,
            name=isolate,
            sequence_key=key,
            subtype=match.group("subtype"),
            isolate_key=isolate,
        )


class UUIDHeader(FastaHeader):
    """Parsing and Validation of exported Fasta names."""

    @classmethod
    def samples(cls) -> list[str]:
        return [
            "5e4399f2-ae07-3fae-e9a8-8e2abb8cbb31/6/H10N1",
        ]

    HEADER_RE: ClassVar[re.Pattern] = re.compile(
        "/".join(  # noqa: FLY002
            [
                RE_UUID,
                RE_SEGMENT_NUM,
                RE_SUBTYPE,
            ]
        )
    )

    @classmethod
    def parse(cls, header: str) -> Self:
        match = cls.do_parse(cls.HEADER_RE, header)
        return cls(
            segment=SegmentType(int(match.group("segment_num"))),
            name=match.group("uuid"),
            sequence_key=match.group("uuid"),
            subtype=match.group("subtype"),
        )
