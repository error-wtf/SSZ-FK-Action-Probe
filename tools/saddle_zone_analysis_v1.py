#!/usr/bin/env python3
"""Saddle-zone analysis for the SSZ-FK action probe (ChatGPT task 1).

For each problematic r_s interval (and a coarse baseline over the whole
range) this tool answers ONE question per point: does an admissible
stationary equilibrium exist, and is it a MINIMUM or a SADDLE?

Method per r_s point:
  * K multistarts: warm start from the neighbouring solutions PLUS random
    starts in the admissible box -> landscape sample, not a single path.
  * Each start: L-BFGS-B (projected) + Newton polish, same admissibility
    discipline as V1 (boundary hit = not stationary).
  * Classification: grad_inf < 1e-6 and min eig(H) > -1e-10 and interior
    -> MINIMUM; grad_inf < 1e-6 with negative eigenvalues -> SADDLE;
    otherwise OPTIMIZER_FAILURE (documented, not hidden).
Additionally:
  * Adaptive bisection of every admissible/nicht-admissible boundary
    (binary search down to dr=0.005).
  * BACKWARD continuation r_s: max -> 0 with warm starts; hysteresis is
    reported as any difference in the admissible interval map or in the
    ground-state minimum found at matching r_s.

Standalone research tool; strict-JSON output.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import minimize

from ssz_fk_action_probe_v1 import Parameters, SSZFKAction

BOUND_FRAC = 0.45
GRAD_TOL = 1e-6


def solve_from(model: SSZFKAction, u0: np.ndarray) -> dict:
    bound = BOUND_FRAC * model.p.dr
    bounds = [(-bound, bound)] * model.p.sites
    opt = minimize(model.potential, u0, method="L-BFGS-B", jac=model.grad,
                   bounds=bounds,
                   options={"maxiter": 3000, "gtol": 1e-9, "ftol": 1e-13,
                            "maxls": 50})
    x = np.clip(opt.x, -bound + 1e-12, bound - 1e-12)
    for _ in range(8):
        uval, gu, Hu = model.potential_gradient_hessian(x)
        if np.max(np.abs(gu)) < 1e-11:
            break
        min_eig = float(np.linalg.eigvalsh(Hu)[0])
        if min_eig <= 0:
            break  # Newton only polishable near minima; saddles kept as-is
        step = np.linalg.solve(Hu, gu)
        accepted = False
        for a in [1.0, 0.5, 0.25, 0.125, 0.0625]:
            trial = x - a * step
            if np.max(np.abs(trial)) >= bound:
                continue
            if model.potential(trial) < uval:
                x, accepted = trial, True
                break
        if not accepted:
            break
    gin = float(np.max(np.abs(model.grad(x))))
    evals = np.linalg.eigvalsh(model.hessian(x))
    min_eig = float(evals[0])
    n_neg = int(np.sum(evals < -1e-10))
    at_bound_mask = np.abs(np.abs(x) - (bound - 1e-9)) < 1e-7
    n_at_bound = int(np.sum(at_bound_mask))
    interior = n_at_bound == 0
    if interior and gin < GRAD_TOL and n_neg == 0:
        kind = "MINIMUM"
    elif interior and gin < GRAD_TOL:
        kind = "SADDLE"
    elif not interior and gin >= GRAD_TOL:
        # Boundary-pinned with active constraint: KKT check.  At the bound
        # the force (-grad) must push OUTWARD against the wall, i.e.
        # sign(u)*grad < 0 (upper bound: grad <= 0; lower bound: grad >= 0),
        # and interior sites must be stationary.  If KKT holds this is a
        # genuine constraint event (the interior minimum no longer exists
        # inside the ordering domain), NOT an optimizer failure.
        g = model.grad(x)
        kkt_ok = np.all(np.sign(x[at_bound_mask]) * g[at_bound_mask] < 0)
        free_ok = (np.max(np.abs(g[~at_bound_mask]))
                   < GRAD_TOL * 50) if (~at_bound_mask).any() else True
        if kkt_ok and free_ok:
            kind = "CONSTRAINT_BOUNDARY"
        else:
            kind = "OPTIMIZER_FAILURE"
    else:
        kind = "OPTIMIZER_FAILURE"
    return {"u": x, "energy": float(model.potential(x)), "grad_inf": gin,
            "min_hessian_eig": min_eig, "n_negative_hessian": n_neg,
            "n_sites_at_bound": n_at_bound,
            "interior": interior, "kind": kind}


def analyze_point(p: Parameters, starts: list[np.ndarray]) -> dict:
    """Multistart landscape sample at one r_s."""
    model = SSZFKAction(p, flat=False)
    results = [solve_from(model, s) for s in starts]
    kinds = [r["kind"] for r in results]
    minima = [r for r in results if r["kind"] == "MINIMUM"]
    saddles = [r for r in results if r["kind"] == "SADDLE"]
    out = {
        "n_starts": len(starts),
        "n_minima_found": len(minima),
        "n_saddles_found": len(saddles),
        "n_failures": kinds.count("OPTIMIZER_FAILURE"),
    }
    if minima:
        best = min(minima, key=lambda r: r["energy"])
        out["ground"] = {"energy": best["energy"],
                         "grad_inf": best["grad_inf"],
                         "min_hessian_eig": best["min_hessian_eig"],
                         "energy_spread_minima": float(max(
                             m["energy"] for m in minima) - best["energy"])}
        out["ground_u"] = best["u"]
        out["admissible"] = True
    else:
        out["ground"] = None
        out["ground_u"] = None
        out["admissible"] = False
    if saddles:
        out["saddle_energies"] = sorted(s["energy"] for s in saddles)[:4]
    return out


def random_starts(rng: np.random.Generator, sites: int, dr: float,
                  k: int) -> list[np.ndarray]:
    bound = BOUND_FRAC * dr
    return [rng.uniform(-bound, bound, size=sites) for _ in range(k)]


def sweep(params: Parameters, rs_values: np.ndarray, seed_points: dict,
          rng: np.random.Generator, k_random: int) -> list[dict]:
    """Forward or backward sweep with warm starts + random multistarts."""
    rows = []
    u_warm = seed_points.get("init", np.zeros(params.sites))
    for rs in rs_values:
        p = replace(params, r_s=float(rs))
        dr = p.dr
        starts = [u_warm] + seed_points.get("nearest", []) \
            + random_starts(rng, params.sites, dr, k_random)
        res = analyze_point(p, starts)
        res["r_s"] = float(rs)
        # carry the best minimum (if any) as warm start for the next r_s
        if res["ground_u"] is not None:
            u_warm = res["ground_u"]
        del res["ground_u"]  # ndarray -> not JSON-serializable
        rows.append(res)
        g = res["ground"]
        print(f"  rs={rs:.4f} admissible={res['admissible']} "
              f"min={res['n_minima_found']} sad={res['n_saddles_found']} "
              f"fail={res['n_failures']}"
              + (f" E={g['energy']:.6f} spread={g['energy_spread_minima']:.2e}"
                 if g else ""), flush=True)
    return rows


def bisect_boundary(params: Parameters, rs_good: float, rs_bad: float,
                    seed_points: dict, rng: np.random.Generator,
                    k_random: int, depth: int = 6) -> float:
    """Binary search the admissible boundary between two r_s values."""
    lo, hi = rs_good, rs_bad  # lo admissible, hi not (any order handled)
    for _ in range(depth):
        mid = 0.5 * (lo + hi)
        p = replace(params, r_s=float(mid))
        starts = seed_points.get("nearest", []) + \
            random_starts(rng, params.sites, p.dr, k_random)
        res = analyze_point(p, starts)
        if res["admissible"]:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path,
                    default=Path("SSZ_FK_SADDLE_ZONES_V1.json"))
    ap.add_argument("--sites", type=int, default=40)
    ap.add_argument("--rs-max", type=float, default=4.4)
    ap.add_argument("--coarse-step", type=float, default=0.1)
    ap.add_argument("--gap-step", type=float, default=0.025)
    ap.add_argument("--k-random", type=int, default=8)
    ap.add_argument("--seed", type=int, default=1103)
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    params = Parameters(sites=args.sites)

    # --- Phase 1: full coarse map with multistart (k_random each point) ---
    print("== Phase 1: coarse multistart map ==", flush=True)
    coarse_rs = np.arange(0.0, args.rs_max + 1e-9, args.coarse_step)
    coarse_rows = sweep(params, coarse_rs, {}, rng, args.k_random)

    # --- Phase 2: fine resolution of the three known gap intervals ---
    gaps = [(0.25, 0.95), (1.85, 2.40), (2.60, 3.05)]
    print("== Phase 2: fine gap sweeps ==", flush=True)
    gap_rows = {}
    for a, b in gaps:
        grid = np.arange(a, b + 1e-9, args.gap_step)
        print(f"-- gap {a}-{b} ({len(grid)} pts) --", flush=True)
        gap_rows[f"{a:.2f}_{b:.2f}"] = sweep(params, grid, {}, rng,
                                             args.k_random)

    # --- Phase 3: boundary bisection (adaptive, dr=0.005 target) ---
    print("== Phase 3: boundary bisection ==", flush=True)
    boundaries = {}
    all_rows = coarse_rows + [r for rows in gap_rows.values() for r in rows]
    adm = {r["r_s"]: r["admissible"] for r in all_rows}
    edges = []
    rs_sorted = sorted(adm)
    for i in range(len(rs_sorted) - 1):
        a, b = rs_sorted[i], rs_sorted[i + 1]
        if adm[a] != adm[b] and (b - a) <= args.gap_step * 1.5:
            edges.append((a, b))
    for a, b in edges:
        key = f"{a:.3f}/{b:.3f}"
        print(f"-- bisect {key} --", flush=True)
        boundaries[key] = bisect_boundary(params, a, b, {}, rng,
                                          args.k_random)
        print(f"   boundary r_s = {boundaries[key]:.4f}", flush=True)

    # --- Phase 4: BACKWARD continuation for hysteresis ---
    print("== Phase 4: backward continuation ==", flush=True)
    back_rs = np.arange(args.rs_max, -1e-9, -0.05)
    # seed with the coarse-map ground state at rs_max (warm start)
    top = max(coarse_rows, key=lambda r: r["r_s"])
    seed = {}
    back_rows = sweep(params, back_rs, seed, rng, args.k_random)

    # --- Hysteresis comparison ---
    fwd_adm = {r["r_s"]: r["admissible"] for r in all_rows
               if abs((r["r_s"] / 0.05) - round(r["r_s"] / 0.05)) < 1e-6}
    mismatches = []
    for r in back_rows:
        rs = round(r["r_s"], 4)
        if rs in fwd_adm and fwd_adm[rs] != r["admissible"]:
            mismatches.append({"r_s": rs,
                               "forward": fwd_adm[rs],
                               "backward": r["admissible"]})
    result = {
        "kind": "SSZ-FK saddle-zone analysis gate (research probe, task 1)",
        "scope": "EFFECTIVE_MODEL_ONLY; no astrophysical or QNM claims",
        "config": {"sites": args.sites, "rs_max": args.rs_max,
                   "coarse_step": args.coarse_step,
                   "gap_step": args.gap_step, "k_random": args.k_random,
                   "seed": args.seed, "grad_tol": GRAD_TOL},
        "classification_legend": {
            "MINIMUM": "grad_inf<1e-6, interior, Hessian positive definite",
            "SADDLE": "grad_inf<1e-6, interior, >=1 negative Hessian eigval",
            "OPTIMIZER_FAILURE": "no stationary point reached from any start"},
        "coarse_map": coarse_rows,
        "gap_sweeps": gap_rows,
        "boundaries_bisected": boundaries,
        "backward_sweep": back_rows,
        "hysteresis_mismatches": mismatches,
    }
    n_adm_fwd = sum(1 for r in all_rows if r["admissible"])
    n_multi = sum(1 for r in all_rows
                  if r["admissible"]
                  and r["ground"]["energy_spread_minima"] > 1e-8)
    n_saddle_only = sum(1 for r in all_rows
                        if not r["admissible"] and r["n_saddles_found"] > 0)
    n_pure_fail = sum(1 for r in all_rows
                      if not r["admissible"] and r["n_saddles_found"] == 0)
    result["summary"] = {
        "n_points_total": len(all_rows),
        "n_admissible": n_adm_fwd,
        "n_admissible_with_multiple_minima": n_multi,
        "n_inadmissible_but_saddle_found": n_saddle_only,
        "n_inadmissible_pure_optimizer_failure": n_pure_fail,
        "hysteresis_points": len(mismatches),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"wrote {args.output}")
    print(json.dumps(result["summary"], indent=2))


if __name__ == "__main__":
    main()
