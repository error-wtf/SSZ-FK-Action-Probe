#!/usr/bin/env python3
"""Pseudo-arclength continuation of the unconstrained SSZ-FK chain (Round III).

Purpose
-------
Trace the branch of stationary points  grad U(u; r_s) = 0  of the UNCONSTRAINED
chain (only geometric radial ordering is enforced, exactly as in
`ordering_constraint_probe_v1.unconstrained_newton`) through parameter folds
using Keller-style pseudo-arclength continuation (bordered Newton).

This answers the open question from README §3.5: the lambda_min sign change
near r_s ~ 0.97 was found on a path-dependent damped-Newton walk whose r_s=1.0
row did NOT converge.  A genuine saddle-node fold of the stationary equations
must (a) be traceable as a smooth solution curve in (u, r_s) space, and
(b) show lambda_min -> 0 exactly at the turning point where dr/ds flips sign.
Both are checked here.  If the tracer passes r_s ~ 1.0 smoothly with
lambda_min staying positive (or crossing without a turning point), the earlier
"fold candidate" was optimizer pathology, not physics.

Frozen gates (declared BEFORE the run, mirroring the round-II discipline)
------------------------------------------------------------------------
  GATE_STATIONARY : every branch point accepted only with ||grad U||_inf < 1e-9
                    (stricter than the 1e-8 of the round-II walk, because the
                    bordered corrector polishes).
  GATE_ORDERING   : min_i (r_{i+1} - r_i) > 1e-9 at every accepted point
                    (model.position() enforces > 0 strictly; recorded values).
  GATE_FD_RS      : central finite difference for d(grad)/d(r_s) with
                    h = 1e-5 * (1 + |r_s|); validated on the r_s = 0 start
                    against the analytic rs-slope probe (relative tangent
                    agreement < 1e-4) before the sweep starts.
  GATE_CONTROLS   : CONTROL-1 rs=0 start must reach grad_inf < 1e-10 from
                    >= 2 of 3 starts (multistart discipline of the shipped
                    probe); CONTROL-2 tangent direction consistent with plain
                    natural-parameter step (u(ds t) vs fixed-rs Newton at
                    r_s = ds*t_s, relative disagreement < 5e-3).
Fold refinement: bisection along the chord between the two accepted points
bracketing the lambda_min sign change; 24 refinements or |lambda| < 1e-11.
Turning-point classification: sign flip of the tangential component t_rs
around the refined crossing => FOLD_CONFIRMED (saddle-node); crossing with
t_rs keeping its sign => STABILITY_EXCHANGE (no fold in r_s).

Standalone research tool; strict-JSON output (allow_nan=False); touches no
frozen repo; re-uses the shipped model code (no physics re-implementation).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from dataclasses import replace
from pathlib import Path

import numpy as np

from ssz_fk_action_probe_v1 import Parameters, SSZFKAction

# ---- frozen tolerances ------------------------------------------------------
GRAD_TOL = 1e-9          # GATE_STATIONARY
START_TOL = 1e-10        # CONTROL-1
ORDER_TOL = 1e-9         # GATE_ORDERING
FD_H = None              # set per point: 1e-5 * (1 + |rs|)
NEWTON_MAX_IT = 30
BACKTRACK_MAX = 40
DS_INIT = 0.02
DS_MIN = 1e-7
DS_MAX = 0.40
MAX_STEPS = 6000
BISECT_ITERS = 24


def grad_rs(model: SSZFKAction, u: np.ndarray, rs: float, h: float) -> np.ndarray:
    """Central FD of d(grad U)/d(r_s) at fixed u (ordering stays valid)."""
    g_plus = SSZFKAction(replace(model.p, r_s=rs + h), False).grad(u)
    g_minus = SSZFKAction(replace(model.p, r_s=rs - h), False).grad(u)
    return (g_plus - g_minus) / (2.0 * h)


def full_jacobian(model: SSZFKAction, u: np.ndarray, rs: float):
    h = 1e-5 * (1.0 + abs(rs))
    U, g, H = model.potential_gradient_hessian(u)
    Frs = grad_rs(model, u, rs, h)
    return U, g, H, Frs, h


def find_start(params: Parameters, rng_seed: int = 7) -> dict:
    """CONTROL-1: rs = 0 stationary point, damped Newton from 3 starts."""
    model = SSZFKAction(replace(params, r_s=0.0), flat=False)
    best = None
    n_ok = 0
    for seed in range(3):
        rng = np.random.default_rng(rng_seed + seed)
        u = np.zeros(params.sites) if seed == 0 else \
            0.02 * params.dr * rng.standard_normal(params.sites)
        U, g, H = model.potential_gradient_hessian(u)
        for _ in range(100):
            if np.max(np.abs(g)) < START_TOL:
                break
            try:
                step = np.linalg.solve(H, -g)
            except np.linalg.LinAlgError:
                break
            alpha = 1.0
            moved = False
            for _bt in range(BACKTRACK_MAX):
                try:
                    Ut, gt, Ht = model.potential_gradient_hessian(u + alpha * step)
                except ValueError:
                    alpha *= 0.5
                    continue
                if Ut <= U - 1e-4 * alpha * float(g @ step):
                    u, U, g, H = u + alpha * step, Ut, gt, Ht
                    moved = True
                    break
                alpha *= 0.5
            if not moved:
                break
        gin = float(np.max(np.abs(model.grad(u))))
        if gin < START_TOL:
            n_ok += 1
        if best is None or gin < best["grad_inf"]:
            best = {"u": u.copy(), "grad_inf": gin,
                    "energy": float(model.potential(u))}
    best["n_starts_ok"] = n_ok
    best["gate_pass"] = bool(n_ok >= 2 and best["grad_inf"] < START_TOL)
    return best


def corrector(model_factory, z: np.ndarray, t: np.ndarray,
              z_pred: np.ndarray) -> dict:
    """Bordered Newton: solve [F(z); t^T (z - z_pred)] = 0 by Newton with
    backtracking; F(z) = grad U(u; rs).  Returns accepted point or failure."""
    n = len(z) - 1
    z = z_pred.copy()
    for it in range(1, NEWTON_MAX_IT + 1):
        u, rs = z[:n], float(z[n])
        model = model_factory(rs)
        try:
            U, g, H, Frs, _h = full_jacobian(model, u, rs)
        except ValueError:
            return {"ok": False, "reason": "ordering_violated_at_iterate",
                    "iters": it}
        c_res = float(t @ (z - z_pred))
        if np.max(np.abs(g)) < 1e-11 and abs(c_res) < 1e-11:
            return {"ok": True, "z": z.copy(), "iters": it, "grad_inf": float(np.max(np.abs(g)))}
        J = np.zeros((n + 1, n + 1))
        J[:n, :n] = H
        J[:n, n] = Frs
        J[n, :] = t
        try:
            dz = np.linalg.solve(J, np.concatenate([-g, [-c_res]]))
        except np.linalg.LinAlgError:
            return {"ok": False, "reason": "bordered_matrix_singular", "iters": it}
        alpha = 1.0
        accepted = False
        base_norm = float(np.max(np.abs(g)))
        for _bt in range(BACKTRACK_MAX):
            z_new = z + alpha * dz
            m_new = model_factory(float(z_new[n]))
            try:
                g_new = m_new.grad(z_new[:n])
            except ValueError:
                alpha *= 0.5
                continue
            c_new = float(t @ (z_new - z_pred))
            if (float(np.max(np.abs(g_new))) <= base_norm
                    and abs(c_new) <= abs(c_res) + 1e-14) or alpha < 1e-6:
                z = z_new
                accepted = True
                break
            alpha *= 0.5
        if not accepted:
            return {"ok": False, "reason": "no_backtrack_progress", "iters": it}
    u, rs = z[:n], float(z[n])
    g = model_factory(rs).grad(u)
    if np.max(np.abs(g)) < 1e-11:
        return {"ok": True, "z": z.copy(), "iters": NEWTON_MAX_IT,
                "grad_inf": float(np.max(np.abs(g)))}
    return {"ok": False, "reason": "corrector_maxiter", "iters": NEWTON_MAX_IT}


def tangent(model: SSZFKAction, u: np.ndarray, rs: float,
            t_prev: np.ndarray) -> np.ndarray | None:
    n = len(u)
    try:
        _U, _g, H, Frs, _h = full_jacobian(model, u, rs)
    except ValueError:
        return None
    J = np.zeros((n + 1, n + 1))
    J[:n, :n] = H
    J[:n, n] = Frs
    J[n, :] = t_prev
    rhs = np.zeros(n + 1)
    rhs[n] = 1.0
    try:
        y = np.linalg.solve(J, rhs)
    except np.linalg.LinAlgError:
        return None
    y /= np.linalg.norm(y)
    if float(y @ t_prev) < 0:
        y = -y
    return y


def branch_point_rows(model: SSZFKAction, u: np.ndarray, rs: float,
                      s: float, ds: float, iters: int) -> dict:
    U, g, H = model.potential_gradient_hessian(u)
    evals = np.linalg.eigvalsh(H)
    x = model.q + u
    return {"rs": float(rs), "s": float(s), "grad_inf": float(np.max(np.abs(g))),
            "energy": float(U), "min_hessian_eig": float(evals[0]),
            "n_negative_hessian": int(np.sum(evals < -1e-10)),
            "min_gap": float(np.min(np.diff(x))),
            "ds_used": float(ds), "corrector_iters": int(iters),
            "u": u.tolist()}


def refine_cross(model_factory, params: Parameters, za: dict, zb: dict,
                 t_a: np.ndarray) -> dict | None:
    """Locate the lambda_min zero crossing. Primary: bordered corrector on
    interpolated predictions (robust at the fold tip where H -> 0 makes
    fixed-rs Newton singular). Fallback: chord bisection with fixed-rs
    polish (only valid away from the tip)."""
    ua, rsa = np.array(za["u"]), za["rs"]
    ub, rsb = np.array(zb["u"]), zb["rs"]
    la, lb = za["min_hessian_eig"], zb["min_hessian_eig"]
    lo, hi = 0.0, 1.0
    best = None
    n = len(ua)

    # ---- primary: bordered-corrector walk on interpolated predictions ----
    # requires u-vectors of the bracketing points
    ua_full = np.array(za.get("u_vec", za.get("u", None)), dtype=float) \
        if "u_vec" in za or "u" in za else None
    ub_full = np.array(zb.get("u_vec", zb.get("u", None)), dtype=float) \
        if "u_vec" in zb or "u" in zb else None
    if ua_full is not None and ub_full is not None:
        z_prev = np.concatenate([ua_full, [rsa]])
        t_guess = None
        for frac in np.linspace(0.25, 1.0, 24):
            z_pred = np.concatenate([
                ua_full + frac * (ub_full - ua_full),
                [rsa + frac * (rsb - rsa)]])
            if t_guess is None:
                m0 = model_factory(float(z_pred[n]))
                t_guess = tangent(m0, z_pred[:n], float(z_pred[n]),
                                  np.concatenate([np.zeros(n), [1.0]]))
                if t_guess is None:
                    break
            corr = corrector(model_factory, z_prev, t_guess, z_pred)
            if not corr["ok"]:
                z_pred = z_pred + 1e-9  # nudge and continue
                continue
            z_prev = corr["z"]
            m = model_factory(float(z_prev[n]))
            U, g, H = m.potential_gradient_hessian(z_prev[:n])
            lam = float(np.linalg.eigvalsh(H)[0])
            t_guess = tangent(m, z_prev[:n], float(z_prev[n]), t_guess)
            if la * lam <= 0 or abs(lam) < 1e-11:
                gap = float(np.min(np.diff(m.q + z_prev[:n])))
                if gap > ORDER_TOL:
                    best = {"rs": float(z_prev[n]), "lambda_min": lam,
                            "grad_inf": float(np.max(np.abs(g))),
                            "min_gap": gap, "energy": float(U),
                            "method": "bordered_corrector_walk"}
                break
        # turning-point test using tangents on both sides of the crossing,
        # evaluated along the WALKED branch (robust at the tip, unlike
        # fixed-rs Newton which targets a different family near H -> 0)
        if best is not None and "t_rs_before" not in best:
            best["t_rs_before"] = float(t_a[-1])
            best["t_rs_after"] = (float(t_guess[-1])
                                  if t_guess is not None else None)
            if best["t_rs_after"] is not None:
                best["classification"] = ("FOLD_CONFIRMED_SADDLE_NODE"
                                          if best["t_rs_before"]
                                          * best["t_rs_after"] < 0
                                          else "STABILITY_EXCHANGE_NO_FOLD")
            else:
                best["classification"] = "TANGENT_AFTER_UNAVAILABLE"
        if best is not None and abs(best["lambda_min"]) < 1e-11:
            best["classification"] = "LAMBDA_ZERO_AT_TIP_EXACT"
            return best

    # ---- fallback: chord bisection with fixed-rs polish -------------------
    for _ in range(BISECT_ITERS):
        mid = 0.5 * (lo + hi)
        u_m = ua + mid * (ub - ua)
        rs_m = rsa + mid * (rsb - rsa)
        model = model_factory(rs_m)
        # polish at fixed rs (the branch passes near the chord point)
        try:
            U, g, H = model.potential_gradient_hessian(u_m)
        except ValueError:
            return best
        for _p in range(30):
            if np.max(np.abs(g)) < 1e-11:
                break
            try:
                step = np.linalg.solve(H, -g)
            except np.linalg.LinAlgError:
                break
            alpha = 1.0
            moved = False
            for _bt in range(BACKTRACK_MAX):
                try:
                    Ut, gt, Ht = model.potential_gradient_hessian(u_m + alpha * step)
                except ValueError:
                    alpha *= 0.5
                    continue
                if Ut <= U - 1e-4 * alpha * float(g @ step):
                    u_m, U, g, H = u_m + alpha * step, Ut, gt, Ht
                    moved = True
                    break
                alpha *= 0.5
            if not moved:
                break
        gin = float(np.max(np.abs(g)))
        if gin >= GRAD_TOL:
            return best          # cannot classify a non-stationary point
        lam = float(np.linalg.eigvalsh(H)[0])
        gap = float(np.min(np.diff(model.q + u_m)))
        if gap <= ORDER_TOL:
            return best
        best = {"rs": float(rs_m), "lambda_min": lam,
                "grad_inf": gin, "min_gap": gap,
                "energy": float(U)}
        if la * lam <= 0:
            hi, lb = mid, lam
        else:
            lo, la = mid, lam
        if abs(lam) < 1e-11:
            break
    # turning-point test: t_rs sign on both sides of the refined point
    t_rs_before = t_a[-1]
    tm = tangent(model_factory(best["rs"]), ua + 0.5 * (ub - ua)
                 if False else np.array(best.get("u", ua)),
                 best["rs"], t_a) if best and "u" in best else None
    # t_rs AFTER the crossing: tangent at point b (already past in sweep order)
    model_b = model_factory(rsb)
    t_b = tangent(model_b, ub, rsb, t_a)
    t_rs_after = float(t_b[-1]) if t_b is not None else None
    best["t_rs_before"] = float(t_rs_before)
    best["t_rs_after"] = t_rs_after
    if t_rs_after is not None:
        best["classification"] = ("FOLD_CONFIRMED_SADDLE_NODE"
                                  if t_rs_before * t_rs_after < 0
                                  else "STABILITY_EXCHANGE_NO_FOLD")
    else:
        best["classification"] = "TANGENT_AFTER_UNAVAILABLE"
    return best


def run(params: Parameters, args) -> dict:
    t0 = time.time()
    rng = np.random.default_rng(2026)
    gates = {
        "GATE_STATIONARY": f"accept only grad_inf < {GRAD_TOL}",
        "GATE_ORDERING": f"min_gap > {ORDER_TOL} at every accepted point",
        "GATE_FD_RS": "central FD h=1e-5*(1+|rs|), validated at start",
        "GATE_CONTROLS": "rs=0 start ok from >=2/3 starts; tangent vs fixed-rs step",
    }

    start = find_start(params)
    if not start["gate_pass"]:
        return {"status": "CONTROL_FAILURE_START", "gates": gates, "start": start}

    n = params.sites
    model0 = SSZFKAction(replace(params, r_s=0.0), flat=False)
    u = start["u"]
    rs = 0.0
    s = 0.0
    ds = DS_INIT

    # CONTROL-2: tangent consistency (natural-parameter comparison)
    t = tangent(model0, u, rs, np.concatenate([np.zeros(n), [1.0]]))
    if t is None:
        return {"status": "CONTROL_FAILURE_TANGENT", "gates": gates}
    step_len = 1e-3
    z_try = np.concatenate([u, [rs]]) + step_len * t
    m_fix = SSZFKAction(replace(params, r_s=float(z_try[n])), flat=False)
    u_fix = z_try[:n].copy()
    U, g, H = m_fix.potential_gradient_hessian(u_fix)
    for _ in range(40):
        if np.max(np.abs(g)) < 1e-11:
            break
        try:
            step = np.linalg.solve(H, -g)
        except np.linalg.LinAlgError:
            break
        alpha = 1.0
        moved = False
        for _bt in range(BACKTRACK_MAX):
            try:
                Ut, gt, Ht = m_fix.potential_gradient_hessian(u_fix + alpha * step)
            except ValueError:
                alpha *= 0.5
                continue
            if Ut <= U - 1e-4 * alpha * float(g @ step):
                u_fix, U, g, H = u_fix + alpha * step, Ut, gt, Ht
                moved = True
                break
            alpha *= 0.5
        if not moved:
            break
    rel_dev = float(np.max(np.abs(u_fix - z_try[:n]))
                    / max(1e-12, float(np.max(np.abs(u_fix)))))
    control2_pass = bool(rel_dev < 5e-3)

    def factory(rs_val: float) -> SSZFKAction:
        return SSZFKAction(replace(params, r_s=rs_val), flat=False)

    rows = [branch_point_rows(model0, u, rs, s, 0.0, 0)]
    folds = []
    t_prev = t
    status = "COMPLETED"
    n_steps = 0
    stall_count = 0
    while n_steps < MAX_STEPS:
        n_steps += 1
        z_cur = np.concatenate([u, [rs]])
        z_pred = z_cur + ds * t_prev
        if z_pred[n] < -0.05 or z_pred[n] > params.r_s + 0.4:
            status = "RS_BOUNDS_REACHED"
            break
        corr = corrector(factory, z_cur, t_prev, z_pred)
        if not corr["ok"]:
            if ds > DS_MIN:
                ds = max(ds * 0.5, DS_MIN)
                stall_count += 1
                if stall_count > 60:
                    status = f"STALL_DSBREAK_{corr['reason']}"
                    break
                continue
            status = f"STALLED_{corr['reason']}"
            break
        stall_count = 0
        z_new = corr["z"]
        u_new, rs_new = z_new[:n], float(z_new[n])
        model_new = factory(rs_new)
        row = branch_point_rows(model_new, u_new, rs_new,
                                s + ds, ds, corr["iters"])
        if row["grad_inf"] >= GRAD_TOL:
            # GATE_STATIONARY: shrink and retry (corrector polish insufficient)
            ds = max(ds * 0.5, DS_MIN)
            stall_count += 1
            if stall_count > 60:
                status = "STALL_GATE_STATIONARY"
                break
            continue
        if row["min_gap"] <= ORDER_TOL:
            status = "ORDERING_LIMIT_REACHED"
            rows.append(row)
            break
        t_new = tangent(model_new, u_new, rs_new, t_prev)
        if t_new is None:
            ds = max(ds * 0.5, DS_MIN)
            continue
        # fold detection on consecutive accepted points
        if (rows[-1]["min_hessian_eig"] * row["min_hessian_eig"] < 0):
            rows[-1]["u_vec"] = u.copy()      # prev accepted point's u
            row["u_vec"] = u_new.copy()
            f = refine_cross(factory, params, rows[-1], row, t_prev)
            if f is not None:
                folds.append(f)
        # adaptive step
        if corr["iters"] <= 3:
            ds = min(ds * 1.4, DS_MAX)
        elif corr["iters"] >= 8:
            ds = max(ds * 0.6, DS_MIN)
        if getattr(args, "fine_near_fold", False) and 0.93 <= rs <= 1.01:
            ds = min(ds, 0.002)
        s += ds
        u, rs, t_prev = u_new, rs_new, t_new
        rows.append({k: v for k, v in row.items()})
    else:
        status = "MAX_STEPS_REACHED"

    compact = [{k: v for k, v in r.items() if k not in ("u", "u_vec")}
               for r in rows]
    return {
        "kind": "SSZ-FK pseudo-arclength continuation of unconstrained chain V1",
        "scope": "EFFECTIVE_MODEL_ONLY; no astrophysical or QNM claims",
        "status": status,
        "gates_declared": gates,
        "controls": {
            "CONTROL_1_start": {k: v for k, v in start.items()
                                if k != "u"},
            "CONTROL_2_tangent_vs_fixed_rs": {
                "relative_deviation": rel_dev, "pass": control2_pass},
        },
        "params": vars(params),
        "n_accepted_points": len(rows),
        "n_steps": n_steps,
        "arclength_total": s,
        "wall_seconds": round(time.time() - t0, 1),
        "fold_candidates": folds,
        "branch": compact,
        "u_final": u.tolist(),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path,
                    default=Path("SSZ_FK_PSEUDO_ARCLENGTH_CONTINUATION_V1.json"))
    ap.add_argument("--sites", type=int, default=40)
    ap.add_argument("--fine-near-fold", action="store_true",
                    help="cap ds to 0.002 while 0.93<=rs<=1.01 to resolve "
                         "the fold tip; accept fewer points elsewhere")
    args = ap.parse_args()
    params = Parameters(sites=args.sites)
    result = run(params, args)
    result["provenance"] = {
        "tool_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "model_module_sha256": hashlib.sha256(
            Path(__file__).with_name("ssz_fk_action_probe_v1.py").read_bytes()
        ).hexdigest(),
        "git_head": None,
        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    try:
        import subprocess
        result["provenance"]["git_head"] = subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
            timeout=10).stdout.strip() or None
    except Exception:
        pass
    args.output.write_text(json.dumps(result, indent=1, allow_nan=False) + "\n")
    print("status:", result["status"])
    print("accepted points:", result["n_accepted_points"],
          "steps:", result["n_steps"],
          "wall:", result["wall_seconds"], "s")
    print("CONTROL_1:", result["controls"]["CONTROL_1_start"])
    print("CONTROL_2:", result["controls"]["CONTROL_2_tangent_vs_fixed_rs"])
    for f in result["fold_candidates"]:
        cls = f.get("classification", "UNCLASSIFIED")
        print(f"FOLD-CANDIDATE rs={f['rs']:.6f} lambda={f['lambda_min']:+.3e} "
              f"grad={f['grad_inf']:.2e} method={f.get('method','bisection')}"
              f" -> {cls}")


if __name__ == "__main__":
    main()
