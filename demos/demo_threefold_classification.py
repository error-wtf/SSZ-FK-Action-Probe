#!/usr/bin/env python3
"""Demo: the threefold classification of r_s points in the SSZ-FK package.

Shows, on the real model, that every r_s in the scan falls into exactly one
of three classes:

  MINIMUM             — interior stationary point, Hessian positive,
                        admissible ground state
  CONSTRAINT_BOUNDARY — NO interior stationary point exists: the optimizer
                        pins onto the box |u| <= 0.45*dr and the KKT sign
                        condition shows the gradient pushes OUTWARD. The
                        equilibrium has left the interior: a fold, not a
                        saddle, not an optimizer failure.
  OPTIMIZER_FAILURE   — (rare, ~0.3%) genuinely unconverged starts.

Run:  /home/error/physics/venv/bin/python demos/demo_threefold_classification.py
Time: ~10 s.
"""
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

PKG = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PKG / "tools"))

from ssz_fk_action_probe_v1 import Parameters, SSZFKAction  # noqa: E402
from curvature_continuation_v1 import warm_equilibrium  # noqa: E402

GRAD_TOL = 1e-6


def classify(rs: float, u_prev: np.ndarray):
    p = Parameters(sites=40, r_s=rs)
    mdl = SSZFKAction(p, flat=False)
    eq = warm_equilibrium(mdl, u_prev)
    x = eq["u"]
    bound = 0.45 * p.dr
    g = mdl.grad(x)
    at_bound = np.abs(np.abs(x) - (bound - 1e-9)) < 1e-7
    kkt_ok = bool(np.all(np.sign(x[at_bound]) * g[at_bound] < 0)) \
        if at_bound.any() else None
    free = ~at_bound
    free_ok = bool(np.all(np.abs(g[free]) < GRAD_TOL * 50))

    if eq["stationary"]:
        return eq, "MINIMUM", at_bound, kkt_ok, free_ok
    if at_bound.any() and (kkt_ok is True) and free_ok:
        return eq, "CONSTRAINT_BOUNDARY", at_bound, kkt_ok, free_ok
    return eq, "OPTIMIZER_FAILURE", at_bound, kkt_ok, free_ok


def main():
    print(__doc__)
    print("\n--- live classification (N=40, warm continuation) ---")
    print(f"{'r_s':>6} {'class':<20} {'#pinned':>8} {'E':>12} "
          f"{'lam_min':>10}")
    u_prev = np.zeros(40)
    for rs in [0.1, 0.5, 0.6, 0.8, 1.0, 2.0, 2.1, 2.5, 2.8, 3.2, 4.0, 4.3]:
        eq, cls, at_bound, kkt_ok, free_ok = classify(rs, u_prev)
        if cls == "MINIMUM":
            u_prev = eq["u"]  # continue only through admissible points
        print(f"{rs:>6.2f} {cls:<20} {int(at_bound.sum()):>8} "
              f"{eq['energy']:>12.6f} {eq['min_hessian_eig']:>10.4f}")

    print("\n--- interpretation ---")
    print("Inside gaps (e.g. rs=0.5-0.85, 1.9-2.3, 2.65-2.95): the ordered")
    print("ground state hits the box boundary with OUTWARD gradient")
    print("(KKT satisfied) -> no interior equilibrium exists there.")
    print("This is an existence boundary (fold), NOT a saddle and NOT an")
    print("optimizer artifact: full artifact SSZ_FK_SADDLE_ZONES_V1.json,")
    print("655/657 CONSTRAINT_BOUNDARY, 0 saddles, 0 hysteresis points.")


if __name__ == "__main__":
    main()
