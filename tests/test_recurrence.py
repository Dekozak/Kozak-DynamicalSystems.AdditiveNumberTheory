import numpy as np
import pytest

from recurrence import (
    CircleRotation,
    FurstenbergCorrespondence,
    SkewProduct,
    ap_free_coloring_exists,
    correlation,
    find_ap,
    has_ap,
    kac_average,
    mono_ap,
    orbit_closure_recurrence,
    square_difference,
    square_return_profile,
    upper_density,
    van_der_waerden_number,
    weyl_sum,
)

GOLDEN = (np.sqrt(5) - 1) / 2


# ---------- Poincare / Kac -------------------------------------------------

def test_kac_lemma_for_rotation():
    rot = CircleRotation(GOLDEN)
    mean_rt, visits = kac_average(rot.step, 0.05, lambda x: x < 0.1, 100_000)
    assert visits > 5000
    assert abs(mean_rt - 10.0) < 0.2  # Kac: E[return time] = 1/mu(A)


def test_skew_product_tracks_n_squared_alpha():
    sk = SkewProduct(GOLDEN)
    orb = sk.orbit((0.0, 0.0), 50)
    n = np.arange(50)
    expected = (n * n * GOLDEN) % 1.0
    assert np.allclose(orb[:, 1], expected, atol=1e-9)


# ---------- multiple recurrence / Szemeredi --------------------------------

def test_correlation_equals_ap_density():
    # S = even numbers in [0,20): c_3(2) counts m with m, m+2, m+4 even
    word = [1 if m % 2 == 0 else 0 for m in range(20)]
    assert correlation(word, 3, 2) == pytest.approx(8 / 20)
    assert correlation(word, 3, 1) == 0.0  # odd gap breaks parity


def test_find_ap_agrees_with_has_ap():
    word = np.zeros(50, dtype=int)
    for m in (3, 10, 17, 24):  # 4-term AP with gap 7
        word[m] = 1
    start, gap = find_ap(word, 4)
    assert (start, gap) == (3, 7)
    assert has_ap(word, 4)
    assert not has_ap(word, 5)


def test_correspondence_dictionary():
    S = range(0, 1000, 3)
    fc = FurstenbergCorrespondence.from_set(S, 1000)
    assert fc.cylinder_measure() == pytest.approx(upper_density(fc.word))
    # multiples of 3: gap 3 gives maximal multiple recurrence
    n, val = fc.best_gap(3, n_max=10)
    assert n == 3 and val > 0.3


# ---------- van der Waerden ------------------------------------------------

def test_w_2_3_equals_9():
    assert ap_free_coloring_exists(8, 2, 3) is not None
    assert ap_free_coloring_exists(9, 2, 3) is None
    assert van_der_waerden_number(2, 3) == 9


def test_mono_ap_detector():
    assert mono_ap([0, 1, 0, 1, 0, 1], 3) == (0, 2)  # positions 0,2,4 colour 0
    witness = ap_free_coloring_exists(8, 2, 3)
    assert mono_ap(witness, 3) is None


def test_orbit_closure_recurrence_reduces_to_mono_ap():
    rng = np.random.default_rng(1)
    col = rng.integers(0, 2, size=200)
    m, n = orbit_closure_recurrence(col, k=3, window=1)
    assert col[m] == col[m + n] == col[m + 2 * n]


# ---------- equidistribution / Sarkozy -------------------------------------

def test_weyl_sums_decay():
    assert abs(weyl_sum(GOLDEN, 100_000, degree=1)) < 0.01
    assert abs(weyl_sum(GOLDEN, 100_000, degree=2)) < 0.05


def test_square_difference_found_in_dense_set():
    S = list(range(0, 500, 2))  # evens: 4 = 2^2 is a square difference
    a, b, n = square_difference(S)
    assert a - b == n * n and a in S and b in S


def test_square_return_profile_positive_somewhere():
    rng = np.random.default_rng(3)
    word = (rng.random(5000) < 0.2).astype(int)
    prof = square_return_profile(word, 20)
    assert prof.max() > 0
