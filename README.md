# Spectral-Rigidity-Ellipse — New Explore

We use the scripts contained herein to **provide numerical evidence for dynamical spectral rigidity among ellipses of various eccentricities**.
The workflow follows the method introduced by **J. De Simoi, V. Kaloshin and Q. Wei** (see the two companion papers listed at the end of this file).

There are **four Python scripts** in the repository root:

1. `ellipse_axes_and_collisions3.py`
2. `paper_method_T_matrix.py`
3. `T_D_and _inverse.py`
4. `truncation_and_visualization.py`

---

## `ellipse_axes_and_collisions3.py`

**Objective –** *Compute the semi-axes of an ellipse for each eccentricity and generate collision-point data for periodic billiard orbits.*

* Fixes the circumference of the ellipse to 1 and scans eccentricities `e` over a user-defined list (default: step 0.01 in \[0, 1)).
* For every `e` it calculates

  * the semi-major axis `a` (via the Complete Elliptic Integral E);
  * the semi-minor axis `b`;
  * Lazutkin collision amplitudes for all periods `q` up to `max_q`.
* Stores results in high-precision text files via **mpmath** (100 decimal digits by default).

**Generates two kinds of files**

* `e_and_semi_axes.txt` – one row per eccentricity (`e`, `a`, `b`)
* `all_periods_<eccentricity>e_col_amplitudes.txt` – collision amplitudes for periods 2 … `max_q`

---

## `paper_method_T_matrix.py`

**Objective –** *Assemble the (infinite) T-matrix truncated to `max_q × max_j` that linearises the isospectral constraints.*

* Loads the collision-point files generated above.
* Evaluates

$$
\mathcal{T}_{q,j} = \sum _{n=0}^{q-1} \frac{\cos(2\pi j x_n^q)}{\mu(x_n^q)}\sin(\phi_n^q)
$$






  (this matches Eq. (3) in De Simoi et al., 2017).
* Results for every eccentricity are cached in a Python dictionary and serialised with **pickle**.

**Generates one file**

* `T_matrices.pkl` – dictionary mapping eccentricity → `numpy.ndarray` of shape (`max_q`, `max_j`)

---

## `T_D_and _inverse.py`

**Objective –** *Construct analytic renormalisation operators **M** and **D** (plus their inverses) and save them for later use.*

* Builds the boolean divisor matrix **M** and the diagonal sinc matrix **D**.
* Computes and stores

  * `matrix_T_D.pkl` – `D · M`
  * `matrix_T_D_inverse.pkl` – `M⁻¹ · D⁻¹`

  using explicit Möbius-inversion formulas.

---

## `truncation_and_visualization.py`

**Objective –** *Analyse the spectral properties of the renormalised operator to infer dynamical spectral rigidity.*

* Forms the composite matrix


  $$\tilde{T} = T_D^{-1} \cdot T$$
  

  loading `T_matrices.pkl` and `matrix_T_D_inverse.pkl`.
* For a user-supplied list of sizes `N` it truncates the leading `N × N` block, computes eigenvalues and norms.
* Saves diagnostic plots to PNG.

**Generates files (per N, per e)**

* `spectrum_<eccentricity>e_<N>.png` – eigenvalues in the complex plane

*Note –* File names adapt automatically to the eccentricity list, the truncation sizes, and the precision settings chosen when the script is run.

---

## References

1. J. De Simoi, V. Kaloshin & Q. Wei — *Dynamical Spectral Rigidity among Z₂-Symmetric Strictly Convex Domains Close to a Circle* (2017)
2. S. Ayub & J. De Simoi — *Numerical Evidence of Dynamical Spectral Rigidity of Ellipses among Smooth Z₂-Symmetric Domains* (2020)
