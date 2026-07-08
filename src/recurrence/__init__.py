"""recurrence: dynamical recurrence as an engine for additive number theory.

The package walks the classical chain of ideas

    Poincare recurrence
        -> Kac's lemma
        -> topological multiple recurrence  (van der Waerden's theorem)
        -> measure-theoretic multiple recurrence + correspondence principle
           (Szemeredi's theorem)
        -> polynomial recurrence  (Furstenberg-Sarkozy square differences)

with executable, testable versions of each step.
"""

from .systems import CircleRotation, DoublingMap, ShiftSystem, SkewProduct, orbit
from .poincare import first_return_times, kac_average, return_times
from .multiple import correlation, correlation_profile, count_aps, find_ap, has_ap
from .correspondence import FurstenbergCorrespondence, upper_density
from .vanderwaerden import (
    ap_free_coloring_exists,
    mono_ap,
    orbit_closure_recurrence,
    van_der_waerden_number,
)
from .equidistribution import (
    discrepancy_star,
    square_difference,
    square_return_profile,
    weyl_sum,
)

__version__ = "0.1.0"
