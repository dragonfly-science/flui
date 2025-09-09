import re

import pytest

from flui.dna import SegmentType, SegmentTypeError
from flui.header import RE_SEGMENT, RE_SUBTYPE, FastaHeader, NCBIHeader, NZHeader


def test_segment_type():
    """Parsing integer and text segments."""
    st = SegmentType.parse("Pb2 ")
    assert st is SegmentType.PB2

    st = SegmentType.parse(" 1 ")
    assert st is SegmentType.PB2

    st = SegmentType.parse(" 3\n   ")
    assert st is SegmentType.PA

    st = SegmentType.parse("8")
    assert st is SegmentType.NS

    with pytest.raises(SegmentTypeError):
        SegmentType.parse("oeoue")

    with pytest.raises(SegmentTypeError):
        SegmentType.parse("10")


def test_segment_parsing():
    re_seg = re.compile(RE_SEGMENT)
    assert re_seg.fullmatch("PB1") is not None
    assert re_seg.fullmatch("PB3") is None


def test_subtype_parsing():
    re_sub = re.compile(RE_SUBTYPE)
    assert re_sub.fullmatch("H1N1") is not None
    assert re_sub.fullmatch("H3N3") is not None
    assert re_sub.fullmatch("H3NX") is not None
    assert re_sub.fullmatch("H3N1,2") is not None
    assert re_sub.fullmatch("MIXED") is not None
    assert re_sub.fullmatch("Some/Name") is None


@pytest.mark.parametrize(
    ("header", "expected"),
    [
        (
            "1|PB2|A/Mallard/NZL/W19_10_79/2019(H3N8)",
            NZHeader(
                segment=SegmentType.PB2,
                name="A/Mallard/NZL/W19_10_79/2019(H3N8)",
                host="Mallard",
                location="NZL",
                sequence_key="W19_10_79/1",
                subtype="H3N8",
                collection_date="2019",
                isolate_key="W19_10_79",
            ),
        ),
        (
            "8|NS|A/Brown_skua/Torgersen_island/o81-b82/2024(H5N1)",
            NZHeader(
                segment=SegmentType.NS,
                name="A/Brown_skua/Torgersen_island/o81-b82/2024(H5N1)",
                host="Brown_skua",
                location="Torgersen_island",
                sequence_key="o81-b82/8",
                subtype="H5N1",
                collection_date="2024",
                isolate_key="o81-b82",
            ),
        ),
    ],
)
def test_nz_header(header: str, expected: NZHeader):
    assert NZHeader.parse(header) == expected


@pytest.mark.parametrize(
    ("header", "expected"),
    [
        (
            "2 | Influenza A virus (A / goose / Guangdong / 1 ) |NC_007358|H5N1",
            NCBIHeader(
                segment=SegmentType.PB1,
                sequence_key="NC_007358",
                name="(A / goose / Guangdong / 1 )",
                subtype="H5N1",
            ),
        ),
    ],
)
def test_ncbi_header(header: str, expected: NCBIHeader):
    """Check that we get rid of whitespace etc."""
    assert NCBIHeader.parse(header) == expected


def test_samples():
    for nm, cls in FastaHeader.registry.items():
        for sample in cls.samples():
            cls.parse(sample)
            assert cls.detect(sample) == nm


def test_detect():
    """Test the detection of the header type."""
    header = "2 | Influenza A virus (A / goose / Guangdong / 1 ) |NC_007358|H5N1"
    assert FastaHeader.detect(header) == "NCBI"
    assert FastaHeader.get_parser(header) == NCBIHeader

    header = "1|PB2|A/Mallard/NZL/W19_10_79/2019(H3N8)"
    assert FastaHeader.detect(header) == "NZ"
    assert FastaHeader.get_parser(header) == NZHeader

    header = "2|PB1|EPI_ISL_131202|A/duck/Alberta/35/1976|CY009609|A/H3N"
    assert FastaHeader.detect(header) == "GISAID"

    header = (
        "2 |Influenza A virus (A/goose/Guangdong/1/1996(H5N1)) polymerase (PB1)"
        "and PB1-F2 protein (PB1-F2) genes, complete cds|NC_007358.1|H5N1"
    )
    assert FastaHeader.detect(header) == "NCBI"
