"""The Furstenberg correspondence principle, made computational.

The correspondence principle is the bridge between combinatorics and
dynamics.  Given a set S of non-negative integers with positive upper
density, form its indicator sequence omega in the full shift {0,1}^N and
let A be the cylinder set { xi : xi_0 = 1 }.  Averaging point masses along
the orbit of omega and passing to a weak-* limit produces a shift-invariant
probability measure mu with

    mu(A) >= upper density of S,

and, crucially,

    mu( A \\cap sigma^{-n} A \\cap ... \\cap sigma^{-(k-1)n} A )
        <= upper density of { m : m, m+n, ..., m+(k-1)n all in S }.

So if Furstenberg's multiple recurrence theorem gives the left-hand side a
positive value for some n, then S contains a k-term arithmetic progression.
Szemeredi's theorem follows.

On a computer we cannot take weak-* limits, but we can compute all the
finite-N objects that the limit is made of, and watch them stabilise.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np

from .multiple import correlation

__all__ = ["FurstenbergCorrespondence", "upper_density"]


def upper_density(word: Iterable[int]) -> float:
    """Density of 1s in a finite word (empirical stand-in for upper density)."""
    w = np.asarray(list(word), dtype=np.int8)
    return float(w.mean()) if w.size else 0.0


@dataclass
class FurstenbergCorrespondence:
    """Finite-N model of the correspondence between a set S and a shift system.

    ``word`` is the indicator sequence of S restricted to [0, N).
    """

    word: np.ndarray

    def __post_init__(self):
        self.word = np.asarray(self.word, dtype=np.int8)

    @classmethod
    def from_set(cls, S: Iterable[int], N: int) -> "FurstenbergCorrespondence":
        w = np.zeros(N, dtype=np.int8)
        for s in S:
            if 0 <= s < N:
                w[s] = 1
        return cls(w)

    # ----- the dictionary, entry by entry ---------------------------------

    def cylinder_measure(self) -> float:
        """Empirical measure of the cylinder [1]: equals the density of S."""
        return upper_density(self.word)

    def multiple_recurrence_measure(self, k: int, n: int) -> float:
        """Empirical mu(A \\cap sigma^{-n}A \\cap ... \\cap sigma^{-(k-1)n}A).

        Positive value  <=>  S contains a k-term AP with gap n.
        """
        return correlation(self.word, k, n)

    def best_gap(self, k: int, n_max: int | None = None) -> tuple[int, float]:
        """Gap n maximising the multiple-recurrence measure, with its value."""
        N = self.word.size
        n_max = n_max or max((N - 1) // max(k - 1, 1), 1)
        best_n, best_v = 0, 0.0
        for n in range(1, n_max + 1):
            v = self.multiple_recurrence_measure(k, n)
            if v > best_v:
                best_n, best_v = n, v
        return best_n, best_v

    def density_along_prefixes(self, num_points: int = 50) -> np.ndarray:
        """Densities over prefixes [0, N_j) -- the sequence whose limsup is
        the upper density, i.e. the quantity the weak-* limit preserves."""
        N = self.word.size
        cuts = np.unique(np.linspace(1, N, num_points, dtype=int))
        csum = np.cumsum(self.word)
        return csum[cuts - 1] / cuts
