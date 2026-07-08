"""Poincare recurrence and Kac's lemma for an irrational circle rotation.

We take T(x) = x + alpha mod 1 with alpha the golden-ratio fractional part,
and A = [0, 0.1).  Poincare says the orbit of any point of A revisits A
forever; Kac says the average waiting time between visits is 1/mu(A) = 10.
"""

import numpy as np

from recurrence import CircleRotation, kac_average, return_times

alpha = (np.sqrt(5) - 1) / 2  # golden rotation, badly approximable
rot = CircleRotation(alpha)
A_len = 0.1
in_A = lambda x: x < A_len

x0 = 0.05  # a point of A
N = 200_000

visits = return_times(rot.step, x0, in_A, N)
mean_rt, n_visits = kac_average(rot.step, x0, in_A, N)

print(f"rotation by alpha = {alpha:.6f}, A = [0, {A_len})")
print(f"visits to A in {N} steps : {n_visits}  (Poincare: infinitely many)")
print(f"first ten return gaps    : {np.diff(visits[:11]).tolist()}")
print(f"empirical mean return    : {mean_rt:.4f}")
print(f"Kac's lemma prediction   : {1 / A_len:.4f}")
