# Recurrence → Arithmetic

This covers what some of my undergraduate research covered supervised by Dr. William Ott

**Recurrence in dynamical systems as an engine for proofs in additive number theory -- made computational.**

One of the most beautiful stories in modern mathematics is that theorems about *orbits of dynamical systems returning near where they started* turn out to be equivalent to theorems about *arithmetic progressions hiding inside sets of integers*. This repository walks that story step by step, with executable, tested Python versions of each link in the chain:

```
Poincaré recurrence  (1890)
   └── Kac's return-time lemma  (1947)
        └── topological multiple recurrence  ⇔  van der Waerden's theorem  (1927 / Furstenberg–Weiss 1978)
             └── measure-theoretic multiple recurrence + correspondence principle
                  ⇔  Szemerédi's theorem  (1975 / Furstenberg 1977)
                       └── polynomial recurrence  ⇔  Furstenberg–Sárközy square differences  (1977–78)
```

None of these theorems can literally be *proved* by a computer running finite experiments, but every object appearing in the proofs -- invariant measures, correlation integrals, orbit closures, return times -- has an exact finite analogue that you can compute, plot, and test. That is what this package does. The point is to make the dictionary between dynamics and combinatorics tangible.

## The dictionary

The Furstenberg correspondence principle translates between the two worlds as follows. Given a set S ⊆ ℕ with positive upper density, let ω ∈ {0,1}^ℕ be its indicator sequence, a point of the *shift space*, and let A be the cylinder set of sequences with a 1 at the origin. Averaging point masses along the orbit of ω under the shift σ and passing to a weak-\* limit produces a shift-invariant probability measure μ. Then:

| Combinatorics (the set S) | Dynamics (the system (X, σ, μ)) |
|---|---|
| upper density of S | μ(A) |
| a return m, m+n ∈ S | a point of A ∩ σ⁻ⁿA |
| a k-term arithmetic progression in S | a point of A ∩ σ⁻ⁿA ∩ ⋯ ∩ σ⁻⁽ᵏ⁻¹⁾ⁿA |
| S − S contains a perfect square | μ(A ∩ σ⁻ⁿ²A) > 0 for some n |

Furstenberg's multiple recurrence theorem says the third row's measure is positive for some n whenever μ(A) > 0 — and via the dictionary this *is* Szemerédi's theorem. The same template with the polynomial n² gives the Furstenberg–Sárközy theorem. Replacing measure by topology (a homeomorphism of a compact space, ε-returns instead of positive measure) gives the Furstenberg–Weiss proof of van der Waerden's theorem, where the compact space is the orbit closure of a colouring inside {1,…,r}^ℤ.

## What's in the package

`recurrence.systems` provides the concrete dynamical systems: the irrational circle rotation x ↦ x+α, the doubling map, the affine skew product (x,y) ↦ (x+α, y+2x+α) on the torus whose orbit of the origin tracks n²α, and shift systems built from integer sets. `recurrence.poincare` computes return times and verifies Kac's lemma (mean return time to A equals 1/μ(A)) along long orbits. `recurrence.multiple` computes the multiple-recurrence correlations c_k(n) = density of m with m, m+n, …, m+(k−1)n all in S — the empirical version of μ(A ∩ T⁻ⁿA ∩ ⋯) — and finds explicit arithmetic progressions. `recurrence.correspondence` packages the finite-N correspondence principle: cylinder measures, best recurrence gaps, and densities along growing prefixes (the raw material of the weak-\* limit). `recurrence.vanderwaerden` contains an exact backtracking decision procedure that certifies small van der Waerden numbers (W(2,3) = 9, W(3,3) = 27) and an orbit-closure search exhibiting the *stronger* recurrence the topological proof actually delivers: shifted colourings agreeing on long windows, not just single positions. `recurrence.equidistribution` verifies the Weyl-sum inputs to polynomial recurrence and finds square differences in explicit positive-density sets.

## Quick start

```bash
git clone https://github.com/Dekozak/Kozak-DynamicalSystems.AdditiveNumberTheory
cd recurrence-arithmetic
pip install -e ".[dev]"
pytest                       # ~10 s: verifies Kac, W(2,3)=9, Weyl decay, ...

python examples/01_poincare_and_kac.py
python examples/02_van_der_waerden.py
python examples/03_szemeredi_correspondence.py
python examples/04_sarkozy_squares.py
```

A taste, from example 3: take the integers with no digit 3 in base 4, thin them by a fair coin, and ask for 4-term progressions. The set has density ≈ 5.6% and no visible structure, yet:

```
set S in [0, 60000): 3343 elements, density 0.0557
multiple recurrence, k = 4:
  best gap n = 2,  empirical mu(A ∩ T⁻ⁿA ∩ T⁻²ⁿA ∩ T⁻³ⁿA) = 0.00290
  explicit 4-term AP in S: [1, 2049, 4097, 6145]  (gap 2048)
```

The explicit progression's gap 2048 = 2·4⁵ is no accident -- the digit condition makes the set nearly invariant under shifts tied to powers of 4, and the search finds this structure on its own. This is a small echo of the structure-versus-randomness dichotomy at the heart of all modern proofs of Szemerédi's theorem.

## Reading list

Furstenberg, *Ergodic behavior of diagonal measures and a theorem of Szemerédi on arithmetic progressions*, J. Analyse Math. 31 (1977). Furstenberg & Weiss, *Topological dynamics and combinatorial number theory*, J. Analyse Math. 34 (1978). Einsiedler & Ward, *Ergodic Theory with a view towards Number Theory*, GTM 259, Springer (2011). Tao, *Poincaré's Legacies* and the blog series on the correspondence principle. McCutcheon, *Elemental Methods in Ergodic Ramsey Theory*, LNM 1722.

## License

MIT
