"""Van der Waerden's theorem: certify W(2,3) = 9 and W(3,3) = 27, then view
the same fact through the Furstenberg-Weiss topological-dynamics lens.
"""

import numpy as np

from recurrence import (
    ap_free_coloring_exists,
    mono_ap,
    orbit_closure_recurrence,
    van_der_waerden_number,
)

# --- exact certification of small van der Waerden numbers -----------------
print("W(2,3) =", van_der_waerden_number(2, 3))   # 9
print("W(3,3) =", van_der_waerden_number(3, 3))   # 27

witness = ap_free_coloring_exists(8, 2, 3)
print("AP-free 2-colouring of {1..8}:", witness)
assert mono_ap(witness, 3) is None
assert ap_free_coloring_exists(9, 2, 3) is None
print("...and none exists for {1..9}: every 2-colouring has a mono 3-AP.\n")

# --- the dynamical viewpoint ----------------------------------------------
# A colouring is a point of the shift space {0,1,2}^Z.  Topological multiple
# recurrence finds shifts sigma^m x, sigma^{m+n} x, ..., agreeing on a long
# window -- much stronger than a single monochromatic AP.
rng = np.random.default_rng(0)
coloring = rng.integers(0, 3, size=5000)

m, n = orbit_closure_recurrence(coloring, k=3, window=3)
print(f"random 3-colouring of length 5000:")
print(f"  shifts by m={m}, m+n, m+2n (n={n}) agree on a window of length 3:")
for j in range(3):
    print(f"  sigma^{m + j * n} x |_[0,3) =", coloring[m + j * n : m + j * n + 3])
print("window length 1 recovers a plain monochromatic 3-term AP.")
