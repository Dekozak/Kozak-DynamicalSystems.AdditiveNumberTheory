"""Van der Waerden's theorem and its topological-dynamics proof, in code.

Van der Waerden (1927): for every r and k there is N = W(r, k) such that
any r-colouring of {1, ..., N} contains a monochromatic k-term arithmetic
progression.

Furstenberg and Weiss (1978) deduced this from *topological multiple
recurrence*: if T is a homeomorphism of a compact metric space X, then for
every k and epsilon > 0 there exist x in X and n >= 1 with

    d(T^{jn} x, x) < epsilon   for j = 1, ..., k-1.

The dictionary: an r-colouring c of the integers is a point of the compact
space {1,...,r}^Z with the shift sigma; two colourings are epsilon-close
when they agree on a long window around the origin; a point x in the orbit
closure of c with sigma^{n}x, sigma^{2n}x, ..., sigma^{(k-1)n}x all close
to x yields integers m, m+n, ..., m+(k-1)n receiving the same colour --
a monochromatic AP.

This module provides (a) an exact backtracking decision procedure for
whether an AP-free colouring of {1..N} exists -- letting you *certify*
small van der Waerden numbers such as W(2,3) = 9 and W(3,3) = 27 -- and
(b) the dynamical translation: locating approximate multiple recurrence
of the shift inside the orbit closure of any concrete colouring.
"""

from __future__ import annotations

from typing import Sequence

import numpy as np

__all__ = [
    "mono_ap",
    "ap_free_coloring_exists",
    "van_der_waerden_number",
    "orbit_closure_recurrence",
]


def mono_ap(coloring: Sequence[int], k: int) -> tuple[int, int] | None:
    """First monochromatic k-term AP in a colouring of {0, ..., N-1}.

    Returns ``(start, gap)`` or ``None``.
    """
    c = list(coloring)
    N = len(c)
    for n in range(1, (N - 1) // (k - 1) + 1):
        for m in range(0, N - (k - 1) * n):
            col = c[m]
            if all(c[m + j * n] == col for j in range(1, k)):
                return m, n
    return None


def _extendable(c: list[int], i: int, k: int) -> bool:
    """Check no monochromatic k-AP *ending* at position i was just created."""
    col = c[i]
    for n in range(1, i // (k - 1) + 1):
        if all(c[i - j * n] == col for j in range(1, k)):
            return False
    return True


def ap_free_coloring_exists(N: int, r: int, k: int) -> list[int] | None:
    """Backtracking search for an r-colouring of {0..N-1} with no
    monochromatic k-term AP.  Returns a witness colouring or ``None``.

    ``None`` for the value N together with a witness for N-1 certifies
    W(r, k) = N (in 1-indexed terms).
    """
    c: list[int] = []

    def rec() -> bool:
        i = len(c)
        if i == N:
            return True
        for col in range(r):
            c.append(col)
            if _extendable(c, i, k) and rec():
                return True
            c.pop()
        return False

    return list(c) if rec() else None


def van_der_waerden_number(r: int, k: int, n_max: int = 200) -> int:
    """Compute W(r, k) exactly by exhaustive search (small cases only).

    W(r, k) is the least N such that every r-colouring of {1..N} has a
    monochromatic k-term AP.  Known small values: W(2,3)=9, W(3,3)=27,
    W(2,4)=35.  Anything larger is out of reach for naive search.
    """
    for N in range(1, n_max + 1):
        if ap_free_coloring_exists(N, r, k) is None:
            return N
    raise RuntimeError(f"W({r},{k}) exceeds search bound {n_max}")


def orbit_closure_recurrence(
    coloring: Sequence[int], k: int, window: int
) -> tuple[int, int] | None:
    """Topological multiple recurrence inside the orbit closure of a colouring.

    Interprets ``coloring`` as a point x of the shift space and searches for
    a base point m and gap n such that the shifted points
    sigma^{m}x, sigma^{m+n}x, ..., sigma^{m+(k-1)n}x pairwise agree on a
    window of length ``window`` -- i.e. they are 2^{-window}-close in the
    usual metric on the shift.  With window = 1 this is exactly a
    monochromatic k-term AP; larger windows exhibit the *stronger* recurrence
    that the dynamical proof actually delivers.

    Returns ``(m, n)`` or ``None`` if the finite sample is too short.
    """
    c = np.asarray(coloring)
    N = c.size
    for n in range(1, N):
        top = N - (k - 1) * n - window
        if top <= 0:
            break
        for m in range(top):
            base = c[m : m + window]
            if all(
                np.array_equal(c[m + j * n : m + j * n + window], base)
                for j in range(1, k)
            ):
                return m, n
    return None
