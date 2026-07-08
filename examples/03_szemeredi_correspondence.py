"""The Furstenberg correspondence principle in action.

Take a positive-density set S, turn it into a point of the shift space,
and read Szemeredi's theorem as multiple recurrence of the cylinder [1]:
some gap n makes  mu(A cap T^-n A cap T^-2n A cap T^-3n A) > 0, which is
literally a 4-term arithmetic progression inside S.

The example set: integers whose base-4 expansion avoids the digit 3, then
thinned by a random 50% coin flip -- density about (3/4)^log-ish, messy,
no obvious structure, yet Szemeredi guarantees long APs.
"""

import numpy as np

from recurrence import FurstenbergCorrespondence, find_ap

N = 60_000
rng = np.random.default_rng(42)


def no_digit_3(m: int) -> bool:
    while m:
        if m % 4 == 3:
            return False
        m //= 4
    return True


S = [m for m in range(N) if no_digit_3(m) and rng.random() < 0.5]
fc = FurstenbergCorrespondence.from_set(S, N)

print(f"set S in [0, {N}): {len(S)} elements, density {fc.cylinder_measure():.4f}")

k = 4
n, val = fc.best_gap(k, n_max=2000)
print(f"\nmultiple recurrence, k = {k}:")
print(f"  best gap n = {n},  empirical mu(A ∩ T^-n A ∩ T^-2n A ∩ T^-3n A) = {val:.5f}")

start, gap = find_ap(fc.word, k)
ap = [start + j * gap for j in range(k)]
print(f"  explicit {k}-term AP in S: {ap}  (gap {gap})")
assert all(a in set(S) for a in ap)

# densities along prefixes -- the sequence whose weak-* limit builds mu
d = fc.density_along_prefixes(8)
print("\ndensities along growing prefixes (the correspondence's raw material):")
print("  ", np.round(d, 4).tolist())
