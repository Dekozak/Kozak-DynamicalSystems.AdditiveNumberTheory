"""Polynomial recurrence: Weyl equidistribution of n^2*alpha and the
Furstenberg-Sarkozy theorem on square differences.

Part 1 checks the dynamical input: for the skew product
T(x,y) = (x+alpha, y+2x+alpha), the orbit of (0,0) has second coordinate
n^2 * alpha, which equidistributes (Weyl) and hence returns near 0.

Part 2 checks the arithmetic output: a dense-ish set with no evident
structure nevertheless contains two elements differing by a perfect square.
"""

import numpy as np

from recurrence import (
    SkewProduct,
    discrepancy_star,
    square_difference,
    square_return_profile,
    weyl_sum,
)

alpha = np.sqrt(2) - 1

# --- Part 1: equidistribution and polynomial recurrence -------------------
for N in (10**3, 10**4, 10**5):
    s1 = abs(weyl_sum(alpha, N, degree=1))
    s2 = abs(weyl_sum(alpha, N, degree=2))
    print(f"N = {N:>6}: |Weyl sum| linear {s1:.5f}   quadratic {s2:.5f}")

sk = SkewProduct(alpha)
orb = sk.orbit((0.0, 0.0), 100_000)
frac = orb[:, 1]  # n^2 * alpha mod 1
print(f"\nstar discrepancy of (n^2 alpha mod 1), N=1e5: {discrepancy_star(frac):.5f}")

near0 = np.flatnonzero(np.minimum(frac, 1 - frac) < 1e-4)[1:6]
print("times n with n^2*alpha within 1e-4 of 0 (polynomial recurrence):", near0.tolist())

# --- Part 2: Sarkozy on a concrete positive-density set -------------------
rng = np.random.default_rng(7)
S = [m for m in range(20_000) if rng.random() < 0.05]  # sparse random, density 5%
a, b, n = square_difference(S)
print(f"\nrandom 5%-density set, {len(S)} elements:")
print(f"  found a = {a}, b = {b} in S with a - b = {n}^2 = {n * n}")

word = np.zeros(20_000, dtype=int)
word[S] = 1
prof = square_return_profile(word, 30)
print("  empirical mu(A ∩ T^-n² A) for n = 1..10:", np.round(prof[:10], 4).tolist())
