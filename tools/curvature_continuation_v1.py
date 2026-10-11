#!/usr/bin/env python3
"""Curvature continuation r_s = 0 -> 4.4 for the action-consistent SSZ-FK probe.

Gate 1 of docs/DERIVATION_AND_FINDINGS.md: track ONE equilibrium branch
continuously in curvature (warm start from the previous r_s solution, no
per-step multistart), and compare each SSZ spectrum against the SAME flat
reference with (a) MAC/Hungarian on individual vectors (V1 diagnostic) and
(b) degeneracy-aware subspace matching via principal angles (the V1
correction).  Records per step: grad_inf, stationarity, lowest Hessian
eigenvalue, negative-Hessian count, alpha fit and residuals.

Standalone research tool; touches nothing in SSZ_FULL_CLOSURE.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import linear_sum_assignment, minimize

from ssz_fk_action_probe_v1 import Parameters, SSZFKAction


def warm_equilibrium(model: SSZFKAction, u0: np.ndarray):
    """Continue the branch: single L-BFGS-B from u0 + Newton polish.

    Same admissibility discipline as the V1 equilibrium(): a constrained
    boundary hit is NOT a stationary point.
    """
    bound = 0.45 * model.p.dr
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
        if np.linalg.eigvalsh(Hu)[0] <= 0:
            break
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
    min_hess = float(np.linalg.eigvalsh(model.hessian(x))[0])
    n_neg = int(np.sum(np.linalg.eigvalsh(model.hessian(x)) < -1e-10))
    return {"u": x, "energy": float(model.potential(x)), "grad_inf": gin,
            "min_hessian_eig": min_hess, "n_negative_hessian": n_neg,
            "stationary": bool(gin < 1e-6 and np.max(np.abs(x)) < bound - 1e-7),
            "at_bound": bool(np.max(np.abs(x)) >= bound - 1e-7)}


def subspace_compare(w_ref, V_ref, mass_ref, w_alt, V_alt, mass_alt,
                     n_top=20, tol=1e-3):
    """Degeneracy-aware frequency comparison.

    Groups near-degenerate frequencies in each spectrum, matches GROUPS via
    principal angles between the mass-weighted eigenvector subspaces, and
    fits alpha + max residual on group-mean frequencies (equal multiplicity
    pairs only, lowest 20).
    """
    # mass-weighted eigenvectors: H v = w2 M v  ->  A = M^{-1/2} H M^{-1/2}
    # eigh already used M; v columns are M-normalized (v^T M v = 1).
    def group_vecs(w, V, mass):
        order = np.argsort(np.where(w > 1e-10, w, np.inf))
        pos = [i for i in order if w[i] > 1e-10]
        groups, cur = [], [pos[0]]
        for i in pos[1:]:
            if abs(np.sqrt(w[i]) - np.sqrt(w[cur[-1]])) <= tol * np.sqrt(w[cur[-1]]):
                cur.append(i)
            else:
                groups.append(cur)
                cur = [i]
        groups.append(cur)
        out = []
        for g in groups:
            # Orthonormal Euclidean basis of the group subspace (same
            # displacement basis for both models -> angles are well defined).
            Qg, _ = np.linalg.qr(V[:, g])
            out.append({"freqs": np.sqrt(np.maximum(w[g], 0)),
                        "mean": float(np.mean(np.sqrt(np.maximum(w[g], 0)))),
                        "Q": Qg})
        return out

    gref, galt = group_vecs(w_ref, V_ref, mass_ref), group_vecs(w_alt, V_alt, mass_alt)
    if not gref or not galt:
        return None
    P = np.array([[np.linalg.svd(a["Q"].T @ b["Q"], compute_uv=False)[0]
                   for b in galt] for a in gref])
    ri, cj = linear_sum_assignment(-P)
    pairs = sorted([(gref[i]["mean"], galt[j]["mean"], float(P[i, j]),
                     len(gref[i]["freqs"]), len(galt[j]["freqs"]))
                    for i, j in zip(ri, cj)])
    same = [p for p in pairs if p[3] == p[4]][:n_top]
    fr = np.array([p[0] for p in same])
    fa = np.array([p[1] for p in same])
    alpha = float(np.dot(fr, fa) / np.dot(fr, fr))
    resid = float(np.max(np.abs(fa - alpha * fr) / (alpha * fr)))
    return {"alpha": alpha, "max_rel_resid": resid, "n_groups": len(same),
            "min_principal_cos": float(np.min([p[2] for p in same])),
            "mean_principal_cos": float(np.mean([p[2] for p in same]))}


def mac_compare(w_ref, V_ref, w_alt, V_alt, n_top=20):
    count = min(n_top, V_ref.shape[1], V_alt.shape[1])
    af = V_ref / np.maximum(np.linalg.norm(V_ref, axis=0), 1e-300)
    bs = V_alt / np.maximum(np.linalg.norm(V_alt, axis=0), 1e-300)
    mac = np.abs(af.T @ bs) ** 2
    row, col = linear_sum_assignment(1 - mac)
    inv = np.empty(V_ref.shape[1], dtype=int)
    inv[row] = col
    ids = np.where(w_ref > 1e-10)[0][:count]
    matched = inv[ids]
    overlap = mac[ids, matched]
    of = np.sqrt(np.maximum(w_ref[ids], 0))
    os_ = np.sqrt(np.maximum(w_alt[matched], 0))
    alpha = float(np.dot(of, os_) / np.dot(of, of))
    resid = float(np.max(np.abs(os_ - alpha * of) / np.maximum(alpha * of, 1e-300)))
    trusted = overlap >= 0.8
    sub = None
    if int(np.sum(trusted)) >= 5:
        at = float(np.dot(of[trusted], os_[trusted]) / np.dot(of[trusted], of[trusted]))
        sub = {"n": int(np.sum(trusted)), "alpha": at,
               "residual": float(np.max(np.abs(os_[trusted] - at * of[trusted])
                                        / np.maximum(at * of[trusted], 1e-300)))}
    return {"alpha": alpha, "residual": resid, "mac_min": float(np.min(overlap)),
            "mac_median": float(np.median(overlap)),
            "n_below_0_5": int(np.sum(overlap < 0.5)), "high_mac": sub}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path,
                    default=Path("SSZ_FK_CURVATURE_CONTINUATION_V1.json"))
    ap.add_argument("--sites", type=int, default=40)
    ap.add_argument("--rs-max", type=float, default=4.4)
    ap.add_argument("--steps", type=int, default=45)
    args = ap.parse_args()

    base_params = Parameters(sites=args.sites)
    steps = np.linspace(0.0, args.rs_max, args.steps)

    flat = SSZFKAction(replace(base_params), flat=True)
    eq0 = warm_equilibrium(flat, np.zeros(base_params.sites))
    w2f, Vf, mf = flat.modes(eq0["u"])

    branch = []
    u_prev = np.zeros(base_params.sites)
    u_prev_ssz = np.zeros(base_params.sites)
    for i, rs in enumerate(steps):
        p = replace(base_params, r_s=float(rs))
        mdl = SSZFKAction(p, flat=False)
        if rs == 0.0:
            eq = eq0
        else:
            eq = warm_equilibrium(mdl, u_prev)
        u_prev = eq["u"]
        record = {"r_s": float(rs), "grad_inf": eq["grad_inf"],
                  "stationary": eq["stationary"], "at_bound": eq["at_bound"],
                  "min_hessian_eig": eq["min_hessian_eig"],
                  "n_negative_hessian": eq["n_negative_hessian"],
                  "energy": eq["energy"]}
        if eq["stationary"] and eq["min_hessian_eig"] > 0:
            w2, V, mass = mdl.modes(eq["u"])
            record["mac"] = mac_compare(w2f, Vf, w2, V)
            sc = subspace_compare(w2f, Vf, mf, w2, V, mass)
            record["subspace"] = sc
            record["rank"] = None
            of = np.sqrt(np.maximum(w2f, 0))
            os_ = np.sqrt(np.maximum(w2, 0))
            ids_f = np.where(w2f > 1e-5)[0][:20]
            ids_s = np.where(w2 > 1e-5)[0][:20]
            if len(ids_f) == len(ids_s):
                a = float(np.dot(of[ids_f], os_[ids_s]) / np.dot(of[ids_f], of[ids_f]))
                record["rank"] = {"alpha": a,
                                  "residual": float(np.max(np.abs(
                                      os_[ids_s] - a * of[ids_f]) / (a * of[ids_f])))}
        branch.append(record)
        r = record.get("rank")
        s = record.get("subspace")
        print(f"[{i+1:02d}/{len(steps)}] rs={rs:.3f} grad={eq['grad_inf']:.2e} "
              f"st={eq['stationary']} negH={eq['n_negative_hessian']} "
              f"rank={('% .4f' % r['residual']) if r else '  n/a '} "
              f"sub={('% .4f' % s['max_rel_resid']) if s else '  n/a '} "
              f"macmin={('% .3f' % record['mac']['mac_min']) if 'mac' in record else 'n/a'}",
              flush=True)

    n_stationary = sum(1 for b in branch if b["stationary"])
    n_unstable = sum(1 for b in branch if b["n_negative_hessian"] > 0)
    result = {"kind": "SSZ-FK curvature continuation gate 1 (research probe)",
              "scope": "EFFECTIVE_MODEL_ONLY; no astrophysical or QNM claims",
              "sites": args.sites, "rs_steps": [float(x) for x in steps],
              "flat_reference_equilibrium": {k: v for k, v in eq0.items() if k != "u"},
              "n_stationary": n_stationary, "n_unstable_steps": n_unstable,
              "branch": branch}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"wrote {args.output}")
    print(f"stationary: {n_stationary}/{len(steps)}, steps with negative "
          f"Hessian eigenvalue: {n_unstable}")


if __name__ == "__main__":
    main()
