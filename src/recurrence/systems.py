"""Concrete dynamical systems used throughout the project.

A (topological / measure-preserving) dynamical system here is simply an
object with a ``step`` method sending a point to its image, together with
helpers to generate orbits.  The systems below are the classical examples
that drive the recurrence-to-arithmetic dictionary:

* ``CircleRotation`` -- the irrational rotation x -> x + alpha (mod 1),
  the simplest ergodic system; underlies Weyl equidistribution and
  three-term recurrence phenomena.
* ``SkewProduct`` -- (x, y) -> (x + alpha, y + 2x + alpha) on the torus,
  whose second coordinate tracks n^2 * alpha.  It is the dynamical engine
  behind Furstenberg/Sarkozy-type theorems ("difference sets of positive
  density sets contain perfect squares").
* ``DoublingMap`` -- x -> 2x (mod 1), a mixing system, for contrast.
* ``ShiftSystem`` -- the left shift on sequences; via the Furstenberg
  correspondence principle it converts subsets of the integers into
  points of a symbolic dynamical system.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from fractions import Fraction
from typing import Callable, Iterator, Sequence

import numpy as np

__all__ = [
    "CircleRotation",
    "SkewProduct",
    "DoublingMap",
    "ShiftSystem",
    "orbit",
]


def orbit(step: Callable, x0, n: int) -> list:
    """Return the finite orbit [x0, T(x0), ..., T^{n-1}(x0)]."""
    xs = [x0]
    for _ in range(n - 1):
        xs.append(step(xs[-1]))
    return xs


@dataclass
class CircleRotation:
    """The rotation T(x) = x + alpha (mod 1) on the circle [0, 1).

    If ``alpha`` is irrational the system is uniquely ergodic with respect
    to Lebesgue measure: orbit averages of any continuous function converge
    to its integral, for every starting point.  This is Weyl's
    equidistribution theorem in dynamical clothing.
    """

    alpha: float

    def step(self, x: float) -> float:
        return (x + self.alpha) % 1.0

    def orbit(self, x0: float, n: int) -> np.ndarray:
        # x_k = x0 + k*alpha mod 1, computed vectorised for accuracy/speed.
        k = np.arange(n)
        return (x0 + k * self.alpha) % 1.0


@dataclass
class SkewProduct:
    """The affine skew product on the 2-torus.

    T(x, y) = (x + alpha, y + 2x + alpha)  (mod 1).

    Starting from (0, 0), the second coordinate after n steps is exactly
    n^2 * alpha (mod 1).  Recurrence of this system therefore says that
    n^2 * alpha returns arbitrarily close to 0 -- the seed of Furstenberg's
    dynamical proof that difference sets of positive-density sets contain
    perfect squares (Sarkozy's theorem).
    """

    alpha: float

    def step(self, p: tuple[float, float]) -> tuple[float, float]:
        x, y = p
        return ((x + self.alpha) % 1.0, (y + 2 * x + self.alpha) % 1.0)

    def orbit(self, p0: tuple[float, float], n: int) -> np.ndarray:
        k = np.arange(n)
        x0, y0 = p0
        xs = (x0 + k * self.alpha) % 1.0
        # closed form: y_k = y0 + k*(2*x0) + k^2*alpha ... derived by summing
        ys = (y0 + 2 * k * x0 + k * k * self.alpha) % 1.0
        return np.column_stack([xs, ys])


@dataclass
class DoublingMap:
    """T(x) = 2x (mod 1): expanding, mixing, positive entropy."""

    def step(self, x: float) -> float:
        return (2.0 * x) % 1.0

    def orbit(self, x0: float, n: int) -> np.ndarray:
        xs = np.empty(n)
        xs[0] = x0
        for i in range(1, n):
            xs[i] = (2.0 * xs[i - 1]) % 1.0
        return xs


@dataclass
class ShiftSystem:
    """The left shift acting on a fixed bi-infinite (here: long finite) word.

    Given the indicator sequence ``word`` of a set S of integers, the point
    omega = word lives in the full shift {0,1}^N and the shift map sigma
    moves the origin one step to the right.  The Furstenberg correspondence
    principle identifies

        density of S            <->  measure of the cylinder [1]
        arithmetic progressions <->  multiple recurrence of that cylinder

    which is how Szemeredi's theorem becomes an ergodic-theoretic statement.
    """

    word: Sequence[int]

    def step(self, i: int) -> int:
        """The shift acts on indices into the word."""
        return i + 1

    def cylinder_indicator(self, i: int) -> int:
        """1 if the shifted point sigma^i(omega) lies in the cylinder [1]."""
        return int(self.word[i]) if 0 <= i < len(self.word) else 0

    def visits(self) -> np.ndarray:
        """Indices i with sigma^i(omega) in the cylinder [1] (i.e. i in S)."""
        w = np.asarray(self.word, dtype=int)
        return np.flatnonzero(w)
