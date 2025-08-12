from collections import defaultdict
from pathlib import Path

import numpy as np
from pydantic import BaseModel, Field, computed_field
from scipy.spatial.distance import jensenshannon

from flui.subtype.counting import (
    MAX_KMER_SIZE,
    MIN_KMER_SIZE,
    int_to_kmer,
    numpy_add_kmers,
)


class KmerSet(BaseModel):
    """A set of related kmers of the same size."""

    size: int = Field(ge=MIN_KMER_SIZE, le=MAX_KMER_SIZE)
    reads: int = 0
    distn: defaultdict[int, int] = Field(
        default_factory=lambda: defaultdict(int), repr=False
    )

    @computed_field
    @property
    def distinct(self) -> int:
        return len(self.distn)

    @computed_field
    @property
    def count(self) -> int:
        return sum(self.distn.values())

    def add_read(self, sequence: bytes, only: set[int] | None = None):
        self.reads += 1
        numpy_add_kmers(self.distn, sequence, self.size, only)

    def add_set(self, other: "KmerSet"):
        if other.size != self.size:
            msg = "Mismatching sizes in kmer sets"
            raise ValueError(msg)

        for k, v in other.distn.items():
            self.distn[k] += v

        self.reads += other.reads

    def jensen_shannon_distance(
        self, other: "KmerSet", only_kmers: set[int] | None = None
    ) -> float:
        """Calculate the Jensen/Shannon Distance.

        We have an option to trim the distribution (use by the index).
        """
        if only_kmers:  # noqa: SIM108
            all_kmers = only_kmers
        else:
            all_kmers = set(self.distn.keys()).union(set(other.distn.keys()))

        prob1 = np.array([self.distn.get(key, 0) for key in all_kmers], dtype=float)
        prob2 = np.array([other.distn.get(key, 0) for key in all_kmers], dtype=float)

        prob1_sum = prob1.sum()
        prob2_sum = prob2.sum()
        if prob1_sum == 0.0 or prob2_sum == 0.0:
            return np.nan

        prob1 /= prob1_sum
        prob2 /= prob2_sum

        # Calculate JS divergence (and convert to python float).
        return float(jensenshannon(prob1, prob2))

    def to_sorted_list(self):
        return sorted(
            [(int_to_kmer(k, self.size), cnt) for (k, cnt) in self.distn.items()]
        )

    def dump(self, out_path: Path):
        with out_path.open(encoding="utf-8", mode="wt") as fd:
            for k, cnt in self.to_sorted_list():
                fd.write(f"{k} {cnt}\n")
