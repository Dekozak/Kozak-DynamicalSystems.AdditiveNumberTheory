"""Poincare recurrence and Kac's return-time lemma, verified numerically.

Poincare (1890): if T preserves a probability measure mu and mu(A) > 0,
then almost every point of A returns to A infinitely often.  Kac (1947)
made this quantitative for ergodic systems: the *expected* first return
time of a point of A is exactly 1 / mu(A).

These two statements are the ancestors of everything else in this project.
Furstenberg's insight was that strengthening "the orbit returns to A" to
"the orbit returns to A along an arithmetic progression of times"
(multiple recurrence) is *equivalent* to Szemeredi's theorem in
combinatorics.
"""

from __future__ import annotations

from typing import Callable, Iterable

import numpy as np

__all__ = ["return_times", "first_return_times", "kac_average"]


def return_times(
    step: Callable, x0, in_A: Callable, n_steps: int
) -> list[int]:
    """Times 0 <= k < n_steps at which the orbit of x0 lies in the set A.

    Parameters
    ----------
    step : the map T
    x0 : starting point
    in_A : indicator function of the set A (returns truthy/falsy)
    n_steps : length of orbit to examine
    """
    times = []
    x = x0
    for k in range(n_steps):
        if in_A(x):
            times.append(k)
        x = step(x)
    return times


def first_return_times(visit_times: Iterable[int]) -> np.ndarray:
    """Gaps between consecutive visits: the sequence of first-return times."""
    t = np.asarray(list(visit_times), dtype=int)
    if t.size < 2:
        return np.array([], dtype=int)
    return np.diff(t)


def kac_average(
    step: Callable, x0, in_A: Callable, n_steps: int
) -> tuple[float, int]:
    """Empirical mean return time to A along one long orbit.

    For an ergodic system Kac's lemma predicts the answer 1 / mu(A).
    Returns ``(mean_return_time, number_of_visits)``.
    """
    visits = return_times(step, x0, in_A, n_steps)
    gaps = first_return_times(visits)
    if gaps.size == 0:
        return float("nan"), len(visits)
    return float(gaps.mean()), len(visits)
