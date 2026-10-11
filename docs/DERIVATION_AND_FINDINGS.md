# Independent action-based SSZ-FK research probe — derivation & findings

**Date:** 2026-10-11. Scope: synthetic, finite effective material lattice; no neutron-star EOS, no real-data claim, no global QNM claim. Standalone experimental variant: original SSZ_FULL_CLOSURE baseline remains frozen.

## 1. Radial action and approximations

For a static radial metric

\[ds^2=-D(r)^2c^2dt^2+s(r)^2dr^2+r^2d\Omega^2,\quad D=(1+r_s/2r)^{-1},\quad s=D^{-1},\]

the coordinate-time action of a radial point particle is

\[L_\mathrm{pp}=-mc^2 D(r)\sqrt{1-\frac{s(r)^2\dot r^2}{D(r)^2c^2}}.\]

For `|s rdot/(D c)| << 1`,

\[L_\mathrm{pp}=-mc^2D(r)+\frac{m s(r)^2}{2D(r)}\dot r^2+O(\dot r^4/c^2).\]

Thus the coefficient of `rdot^2/2` in **coordinate time** is `m*s²/D`; the old `m/D²` is a separately testable effective ansatz, not this point-particle limit.

A spring or pinning energy `U_local` measured in a static local frame contributes `D*U_local` to the coordinate-time potential, subject to the chosen static slicing, matter coupling and neglect of nonlocal retardation. A bond's lapse is approximated by `D((r_i+r_{i+1})/2)`; this is an explicit discretization prescription, **not a unique consequence of relativity**.

For the selected g1 metric,

\[\ell(r)=\int_R^r s(\rho)\,d\rho = (r-R)+\frac{r_s}{2}\log\frac rR,\quad \ell'(r)=s(r).\]

We define `q_i=R+i*dr`, `r_i=q_i+u_i`, the same radial coordinate for metric and material positions. The 1D effective potential is

\[
\begin{aligned}
U &= \sum_i mc^2[D(r_i)-D(q_i)-D'(q_i)u_i]\\
&+\sum_i\frac\kappa2 D((r_i+r_{i+1})/2)\,[\ell(r_{i+1})-\ell(r_i)-\ell^0_i]^2\\
&+\sum_i D(r_i)V_0[1-\cos(Q\ell(r_i))].
\end{aligned}
\]

- The subtracted linear gravity term is an **external support** force at q, mimicking static balance with other matter; it is not derived from the SSZ field equations or a neutron-star stress tensor. Its residual `D''` gives a tidal potential and may destabilize the unpinned translational mode.
- Default natural bond length `ell^0_i=ell(q_{i+1})-ell(q_i)` is locally unstrained at the reference. `constant` is an alternative, potentially prestressed prescription. A constrained energy optimizer may appear successful without finding an unconstrained stationary point: such runs are marked **inadmissible**.
- `Q=2*pi/[dr*(1+delta)]`, delta=(sqrt(5)-2)/4. This remains an assumed FK period; it is **not derived** from a microscopic neutron-star crystal calculation.

The linearized modes satisfy `H v_n = omega_n² diag(mu_i) v_n`, with `mu_i = m s(r_i)^2 / D(r_i)` evaluated at the equilibrium configuration. Normalization `v_i.T @ M @ v_j = delta_ij` implies the source-to-detector resolvent

\[R_{db}(z)= d^T(H-zM)^{-1}b=\sum_n\frac{(d^Tv_n)(v_n^Tb)}{\omega_n^2-z}\]

and sum rule `sum_n |v_n.T @ b|² = b.T @ M^{-1} @ b` for real source `b`. This is what the implementation tests. This finite-Dirichlet-like material response is **not** the outgoing SSZ 3-DOF QNM resolvent.

## 2. Tested phase variants

1. `proper_integral`: `phi=Q ell(r)`; local coordinate wavenumber `Q s(r)`.
2. `literal_sr`: `phi=Q[s(r)r-s(R)R]`. Since `s(r)r=r+r_s/2`, **phi=Q(r-R), derivative Q constant** in the exterior branch.
3. `frozen_local`: `phi=Q s(q_i)(r_i-R)`, a site-dependent effective phase (akin to the old two-scale assumption after choosing a mapping). It is not an integrated proper phase.

There is no mathematical basis for treating all three as physically identical. They are competing prescriptions to be tested separately.

## 3. Numerical pilot results (independent implementation)

Default `R=12`, `r_s=4.4`, `dr=0.05`, `kappa=25`, `m=1`, `V0=0.0003`, `c=1`, natural lengths at reference, with 4 deterministic multistarts and positive-Hessian Newton polishing. All units are **model units**, with no mapping to actual astrophysical Hz. `alpha` is fitted separately in each case. The comparison selects the lowest 20 positive model frequencies and also uses global Hungarian assignment of Euclidean-normalized displacement eigenvectors (mode assurance criterion, MAC).

| Sites | All positive modes: rank-fit max fractional residual | Highest-overlap subset MAC>=0.8: max fractional residual | Min assigned MAC | Interpretation |
|---|---:|---:|---:|---|
| 20 | 0.5785 | 0.4968 (15 modes) | 0.429 | Mode identity ambiguous |
| 40 | 1.1710 | 0.02843 (12 modes) | 0.411 | Low-MAC sector dominates extremes |
| 60 | 1.6381 | 0.01632 (8 modes) | 0.213 | Low-MAC sector strongly affects result |

**Critical finding:** High headline residuals do **not** by themselves establish an invariant spectral effect. At 40/60 sites, restricting to high-MAC modes cuts the max fractional discrepancy to a few percent or less. The low-MAC sector may contain genuine physical rearrangements, different local minima, or identification ambiguity. Neither option has been proved. The MAC threshold is diagnostic, not a preregistered significance threshold; its varying sample size (15,12,8) prohibits a continuity claim.

The 20/40/60 runs change the *number and span* of physical sites, not a fixed-domain discretization, so **these are finite-size sensitivity checks, not continuum convergence**. Additionally, at least some unstable/constraint-bound controls must be excluded from physical comparison.

### Ablations (40 sites)

| Variant | Rank-fit max residual | Min assigned MAC | Qualification |
|---|---:|---:|---|
| proper integral (main) | 1.1710 | 0.411 | ambiguous low-overlap modes |
| literal `s(r)r` phase | 0.32035 | 0.841 | unambiguous at the chosen threshold |
| frozen local site phase | 0.95180 | 0.423 | ambiguous low-overlap modes |
| legacy inertia `m/D²`, same action-potential | 1.17089 | 0.412 | *not* derived kinetic term |
| V0=0, tidal gravity retained | 0.003496 | 0.973 | **negative generalized eigenvalue**; free chain unstable due tidal/support model |
| V0=0, no tidal gravity | 1.11e-5 | 0.973 | admissible; shows near-uniform frequency scaling |
| constant proper natural bond length | 0.81868 | 0.248 | inadmissible: constrained boundary & unstable Hessian |

These are not equal-parameter comparisons to the previous original SSZ FK Run 11; no claim that the new numbers supersede its 0.5642467 residual. They suggest the new signal's attribution is sensitive to the pinning phase and to identity of soft modes. A paper-ready conclusion **cannot** be drawn yet.

## 4. Tests actually run

`pytest -q tests/`: **12 tests passed**, covering analytical `ell'`, three phase derivatives, flat metric limit, proper radial point particle quadratic expansion, six combinations of potential gradient/Hessian compared to independent finite differences, mass-weighted eigensystem and source-to-detector spectral decomposition, and Hungarian assignment with permuted modes. (The test suite consists of 12 collected cases; independent execution, not formal Lean proof.)

## 5. Next scientific gates (not yet passed)

1. **Stationary-state phase continuation:** track a *single* physical equilibrium branch continuously in curvature `r_s` (not separate multistart minima). Report energy barriers, number of negative Hessian eigenvalues, mode branch stability and crossings.
2. **Subspace-aware matching:** degenerate/near-degenerate eigenvalue groups must be matched using principal angles/projectors rather than one-to-one eigenvectors; compare coherent source residues as well as frequencies.
3. **True fixed-domain refinement:** specify continuum elastic/line-density scaling of `m`, `kappa`, `V0` and `dr` first; then fix radial length and refine `n`. Do not call changing the number of physical crystals a grid-convergence theorem.
4. **Crust source physics:** replace assumed pinning with a justified local microphysical elastic / magnetic model and sources `C_J(omega)`; compare standard-GR background using the **same** material action, not just `flat` versus SSZ.
5. **Operator bridge and Lean:** hold the frozen 3DOF SSZ operator fixed until B6 is closed; separately formalize geometric identities, positive-definite mass-matrix normalization, and quadratic response identities in Lean 4. No checked Lean proofs are included in this packet.

## 6. Scientific context

- Upstream canonical branch: `error-wtf/SSZ_FULL_CLOSURE`, `spectroscopy-real-data-20261004`, checked against `ef4dfd3a9efa98282ca237211afb40d303885e7f`.
- Upstream reported FK Run 11: `alpha_best=0.91919070778`, `max_rel_resid_after_alpha_fit=0.5642466802`. These belong to the prior effective two-coordinate coupling, not this probe.
- Existing `PEAKS_NOT_ROBUST` and globally held B6/QNM status are unaffected.
- This research packet is **standalone**; it has not been committed to your protected GitHub branch. It is intended for review before integrating into a separate research branch.
