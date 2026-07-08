"""Multiple recurrence and its equivalence with Szemeredi's theorem.

Furstenberg's multiple recurrence theorem (1977):

    If T preserves a probability measure mu and mu(A) > 0, then for
    every k >= 1 there exists n >= 1 with

        mu( A  \\cap  T^{-n} A  \\cap  T^{-2n} A  \\cap ... \\cap  T^{-(k-1)n} A ) > 0.

Through the correspondence principle (see ``correspondence.py``) this is
*equivalent* to Szemeredi's theorem: every subset of the integers with
positive upper density contains arbitrarily long arithmetic progressions.

This module works with the combinatorial shadow of the measure-theoretic
statement.  For a finite 0/1 word ``w`` (the indicator of a set S inside
an interval [0, N)), the *k-th correlation at gap n*

    c_k(n) = (1/N) * #{ m : w[m] = w[m+n] = ... = w[m+(k-1)n] = 1 }

is precisely the empirical version of mu(A \\cap T^{-n}A \\cap ... ).
Multiple recurrence says: some n makes c_k(n) positive.  Szemeredi says:
S contains a k-term arithmetic progression.  Same statement.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np

__all__ = [
    "correlation",
    "correlation_profile",
    "find_ap",
    "count_aps",
    "has_ap",
]


def correlation(word: Iterable[int], k: int, n: int) -> float:
    """Empirical multiple-recurrence correlation c_k(n) of a 0/1 word."""
    w = np.asarray(list(word), dtype=np.int8)
    N = w.size
    if k < 1 or n < 1 or (k - 1) * n >= N:
        return 0.0
    prod = w[: N - (k - 1) * n].copy()
    for j in range(1, k):
        prod &= w[j * n : N - (k - 1 - j) * n]
    return float(prod.sum()) / N


def correlation_profile(word: Iterable[int], k: int, n_max: int) -> np.ndarray:
    """Array [c_k(1), ..., c_k(n_max)] of correlations at every gap."""
    w = list(word)
    return np.array([correlation(w, k, n) for n in range(1, n_max + 1)])


def find_ap(word: Iterable[int], k: int) -> tuple[int, int] | None:
    """First k-term arithmetic progression inside the set S = {m : w[m]=1}.

    Returns ``(start, gap)`` or ``None``.  This is the combinatorial
    counterpart of exhibiting a point of A that returns to A along the
    progression of times n, 2n, ..., (k-1)n.
    """
    w = np.asarray(list(word), dtype=np.int8)
    N = w.size
    support = np.flatnonzero(w)
    sset = set(support.tolist())
    for m in support:
        max_n = (N - 1 - m) // (k - 1) if k > 1 else 0
        for n in range(1, max_n + 1):
            if all((m + j * n) in sset for j in range(1, k)):
                return int(m), int(n)
    return None


def has_ap(word: Iterable[int], k: int) -> bool:
    """True iff the set encoded by ``word`` contains a k-term AP."""
    return find_ap(word, k) is not None


def count_aps(word: Iterable[int], k: int) -> int:
    """Number of k-term APs in the set (counted by (start, gap) pairs)."""
    w = np.asarray(list(word), dtype=np.int8)
    N = w.size
    total = 0
    for n in range(1, (N - 1) // max(k - 1, 1) + 1):
        prod = w[: N - (k - 1) * n].copy()
        for j in range(1, k):
            prod &= w[j * n : N - (k - 1 - j) * n]
        total += int(prod.sum())
    return total
