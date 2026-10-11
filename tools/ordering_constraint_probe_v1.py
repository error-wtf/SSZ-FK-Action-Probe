#!/usr/bin/env python3
"""Ordering-constraint probe V2 (ChatGPT round-II correction).

Replaces the per-site box |u_i| <= 0.45*dr (sufficient but NOT necessary
for ordering) with the true geometric ordering constraint

    r_{i+1} - r_i >= eps      (i.e.  u_{i+1} - u_i >= eps - dr)

parametrised by gap variables  v = (u_0, g_1, ..., g_{n-1})  with
u_i = u_0 + sum_{j<=i} g_j, so the constraint set becomes a simple box in
v-space and L-BFGS-B applies directly (gradient via chain rule, exact).

Three coordinated tests per r_s point:
  1. GAP-CONSTRAINED optimisation (multistart) with KKT classification:
     MINIMUM (all gaps slack, Hessian PD) / ORDERING_CONTACT (active gap,
     KKT sign OK, free sites stationary) / SADDLE / OPTIMIZER_FAILURE.
  2. UNCONSTRAINED Newton continuation (no constraints at all except
     radial ordering of the evaluation points) with Hessian-spectrum
     monitoring: a lambda_min zero-crossing locates a genuine fold
     candidate of the stationary equations.
  3. COMPARISON against the shipped box-constrained result
     (SSZ_FK_SADDLE_ZONES_V1.json) on the five problematic intervals.

eps is varied (eps_frac * dr in {0.05, 0.2, 0.5}) to test sensitivity.
A regime change that survives the geometric constraint (and matches the
unconstrained lambda_min crossing) earns a physical fold interpretation;
if it does not, the box was the binding structure (modelling artefact).

Standalone research tool; strict-JSON output; touches no frozen repo.
"""
from __future__ import annotations

import argparse
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from ssz_fk_action_probe_v1 import Parameters, SSZFKAction

GRAD_TOL = 1e-6
Hess_TOL = -1e-10
V0_BOUND = 2.0          # bound on the collective shift u_0 (pinning anchors it)
GAP_UPPER = 0.5         # generous upper gap bound (10*dr); ordering needs no upper bound physically


class GapParametrization:
    """v-space wrapper: v = (u_0, g_1..g_{n-1}), u_i = u_0 + cumsum(g)_i."""

    def __init__(self, model: SSZFKAction, eps: float):
        self.model = model
        self.eps = eps
        self.n = model.p.sites
        self.dr = model.p.dr

    def u_of_v(self, v: np.ndarray) -> np.ndarray:
        return v[0] + np.concatenate([[0.0], np.cumsum(v[1:])])

    def potential(self, v):
        return self.model.potential(self.u_of_v(v))

    def grad(self, v):
        gu = self.model.grad(self.u_of_v(v))
        dv0 = float(gu.sum())
        dg = np.array([gu[k + 1:].sum() for k in range(self.n - 1)])
        return np.concatenate([[dv0], dg])

    def hessian_v(self, v):
        u = self.u_of_v(v)
        Hu = self.model.hessian(u)
        J = np.tril(np.ones((self.n, self.n)))
        # v0 column is all-ones; gap columns are step functions
        J = np.hstack([np.ones((self.n, 1)), J[:, 1:]])
        return J.T @ Hu @ J

    def bounds(self):
        lo = self.eps - self.dr
        return ([(-V0_BOUND, V0_BOUND)]
                + [(lo + 1e-12, GAP_UPPER)] * (self.n - 1))

    def solve(self, v_init: np.ndarray) -> dict:
        opt = minimize(self.potential, v_init, method="L-BFGS-B",
                       jac=self.grad, bounds=self.bounds(),
                       options={"maxiter": 4000, "gtol": 1e-10,
                                "ftol": 1e-15, "maxls": 60})
        v = np.clip(opt.x, [b[0] + 1e-12 for b in self.bounds()],
                    [b[1] - 1e-12 for b in self.bounds()])
        u = self.u_of_v(v)
        gv = self.grad(v)
        U, gu, Hu = self.model.potential_gradient_hessian(u)
        evals = np.linalg.eigvalsh(Hu)
        min_eig = float(evals[0])
        n_neg = int(np.sum(evals < Hess_TOL))
        dr = self.dr
        gap_lo = self.eps - dr
        active = np.where(v[1:] <= gap_lo + 1e-8)[0]      # gap indices at lower bound
        u0_at_bound = abs(v[0]) >= V0_BOUND - 1e-7
        gin = float(np.max(np.abs(gv)))
        free_ok = True
        kkt_ok = True
        if len(active):
            # KKT at a lower bound: dU/dg_k >= -tol (pushing outward = down
            # would lower U otherwise).  Free variables near-stationary.
            kkt_ok = bool(np.all(gv[1 + active] >= -1e-7))
            free_mask = np.ones(len(v), bool)
            free_mask[1 + active] = False
            free_mask[0] = not u0_at_bound
            if free_mask.any():
                free_ok = bool(np.max(np.abs(gv[free_mask])) < GRAD_TOL * 50)
        if len(active) == 0 and not u0_at_bound and gin < GRAD_TOL:
            kind = "MINIMUM" if n_neg == 0 else "SADDLE"
        elif gin < GRAD_TOL * 1.0 and kkt_ok and free_ok:
            kind = "ORDERING_CONTACT" if (active or u0_at_bound) else "MINIMUM"
        elif (active or u0_at_bound) and kkt_ok and free_ok:
            kind = "ORDERING_CONTACT"
        else:
            kind = "OPTIMIZER_FAILURE"
        return {"u": u, "v": v, "energy": float(U), "grad_inf": gin,
                "min_hessian_eig": min_eig, "n_negative_hessian": n_neg,
                "n_active_gaps": int(len(active)),
                "active_gap_sites": [int(a) for a in active],
                "u0_at_bound": bool(u0_at_bound), "kind": kind}


def u_to_v(u: np.ndarray) -> np.ndarray:
    return np.concatenate([[u[0]], np.diff(u)])


def random_starts_v(rng, eps, dr, n, k):
    lo = eps - dr + 1e-9
    out = []
    for _ in range(k):
        v0 = rng.uniform(-0.45 * dr, 0.45 * dr)
        gaps = rng.uniform(max(lo, -0.2 * dr), 0.45 * dr, size=n - 1)
        out.append(np.concatenate([[v0], gaps]))
    return out


def analyze_point(params: Parameters, eps: float, warm_u, rng, k_random):
    model = SSZFKAction(params, flat=False)
    gp = GapParametrization(model, eps)
    starts = [u_to_v(warm_u)] + random_starts_v(rng, eps, params.dr,
                                                params.sites, k_random)
    results = [gp.solve(s) for s in starts]
    minima = [r for r in results if r["kind"] == "MINIMUM"]
    contacts = [r for r in results if r["kind"] == "ORDERING_CONTACT"]
    kinds = [r["kind"] for r in results]
    out = {"n_starts": len(starts),
           "n_minima": len(minima),
           "n_ordering_contacts": len(contacts),
           "n_failures": kinds.count("OPTIMIZER_FAILURE"),
           "n_saddles": kinds.count("SADDLE")}
    if minima:
        best = min(minima, key=lambda r: r["energy"])
        out.update({"admissible": True,
                    "energy": best["energy"],
                    "grad_inf": best["grad_inf"],
                    "min_hessian_eig": best["min_hessian_eig"],
                    "energy_spread": float(max(m["energy"] for m in minima)
                                           - best["energy"]),
                    "u": best["u"].tolist()})
    elif contacts:
        best = min(contacts, key=lambda r: r["energy"])
        out.update({"admissible": False,
                    "energy": best["energy"],
                    "grad_inf": best["grad_inf"],
                    "min_hessian_eig": best["min_hessian_eig"],
                    "n_active_gaps": best["n_active_gaps"],
                    "active_gap_sites": best["active_gap_sites"],
                    "u": best["u"].tolist()})
    else:
        out.update({"admissible": False, "energy": None,
                    "min_hessian_eig": None, "u": None})
    out["warm_next"] = (np.array(out["u"]) if out["u"] is not None
                        else warm_u)
    return out


def unconstrained_newton(model: SSZFKAction, u0: np.ndarray,
                         max_iter: int = 80) -> dict:
    """Damped Newton on grad U = 0 with NO constraints; ordering of the
    evaluation points is enforced only by step damping (backtracking)."""
    x = u0.copy()
    U, g, H = model.potential_gradient_hessian(x)
    for _ in range(max_iter):
        if np.max(np.abs(g)) < 1e-10:
            break
        try:
            step = np.linalg.solve(H, -g)
        except np.linalg.LinAlgError:
            return {"converged": False, "reason": "singular Hessian",
                    "u": x.tolist()}
        alpha = 1.0
        ok = False
        for _ in range(30):
            trial = x + alpha * step
            try:
                Ut = model.potential_gradient_hessian(trial)[0]
            except ValueError:
                alpha *= 0.5
                continue
            if Ut <= U - 1e-4 * alpha * float(g @ step):
                x, U, ok = trial, Ut, True
                break
            alpha *= 0.5
        if not ok:
            break
        _, g, H = model.potential_gradient_hessian(x)
    U, g, H = model.potential_gradient_hessian(x)
    evals = np.linalg.eigvalsh(H)
    gin = float(np.max(np.abs(g)))
    return {"converged": bool(gin < 1e-8),
            "grad_inf": gin, "energy": float(U),
            "min_hessian_eig": float(evals[0]),
            "n_negative_hessian": int(np.sum(evals < Hess_TOL)),
            "min_gap": float(np.min(np.diff(model.q + x))),
            "u": x.tolist()}


def bisect_edge(params, eps, rs_good, rs_bad, warm_u, rng, k, depth=7):
    lo, hi = rs_good, rs_bad
    w = warm_u
    for _ in range(depth):
        mid = 0.5 * (lo + hi)
        res = analyze_point(replace(params, r_s=float(mid)), eps, w, rng, k)
        w = res["warm_next"]
        if res["admissible"]:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output", type=Path,
                    default=Path("SSZ_FK_ORDERING_CONSTRAINT_V1.json"))
    ap.add_argument("--sites", type=int, default=40)
    ap.add_argument("--k-random", type=int, default=4)
    ap.add_argument("--seed", type=int, default=1104)
    ap.add_argument("--eps-fracs", type=float, nargs="+",
                    default=[0.05, 0.2, 0.5])
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)
    params = Parameters(sites=args.sites)
    dr = params.dr

    # The five problematic intervals from the box-constrained study.
    intervals = {"gap1": (0.2795, 0.8998), "gap2": (1.8580, 2.3498),
                 "gap3": (2.6057, 2.9502), "endzone_a": (3.35, 3.45),
                 "endzone_b": (4.05, 4.10)}

    box_art = Path("SSZ_FK_SADDLE_ZONES_V1.json")
    box_summary = None
    if box_art.exists():
        bz = json.loads(box_art.read_text())
        box_summary = bz.get("summary")

    # ---- unconstrained Newton continuation (fold detection) ----
    print("== unconstrained Newton continuation (lambda_min monitoring) ==",
          flush=True)
    grid = np.arange(0.0, 4.4 + 1e-9, 0.05)
    unc_rows = []
    u_prev = np.zeros(params.sites)
    flat0 = SSZFKAction(replace(params, r_s=0.0), flat=False)
    warm = unconstrained_newton(flat0, u_prev)
    prev_eig = None
    crossings = []
    for rs in grid:
        mdl = SSZFKAction(replace(params, r_s=float(rs)), flat=False)
        res = unconstrained_newton(mdl, np.array(warm["u"]))
        res["r_s"] = float(rs)
        unc_rows.append({k: v for k, v in res.items() if k != "u"})
        lam = res["min_hessian_eig"]
        if prev_eig is not None and res["converged"]:
            if prev_eig[1] > 0 >= lam or (prev_eig[1] < 0 <= lam):
                a, b = prev_eig[0], float(rs)
                for _ in range(8):
                    m = 0.5 * (a + b)
                    rm = unconstrained_newton(
                        SSZFKAction(replace(params, r_s=float(m)), flat=False),
                        np.array(warm["u"]))
                    if not rm["converged"]:
                        break
                    if (rm["min_hessian_eig"] > 0) == (lam > 0):
                        b = m
                    else:
                        a = m
                crossings.append({"lambda_zero_crossing_rs": 0.5 * (a + b),
                                  "direction": "stability_loss" if lam < 0
                                  else "stability_regain"})
        if res["converged"]:
            warm = res
            prev_eig = (float(rs), lam)
        print(f"  rs={rs:.2f} conv={res['converged']} lam={lam:.3e} "
              f"min_gap={res['min_gap']:.4f}", flush=True)

    # ---- gap-constrained study per eps ----
    per_eps = {}
    for eps_frac in args.eps_fracs:
        eps = eps_frac * dr
        print(f"\n== eps = {eps_frac}*dr = {eps:.4f} ==", flush=True)
        iv_out = {}
        for name, (a, b) in intervals.items():
            grid_iv = np.arange(a - 0.05, b + 0.0501, 0.05)
            grid_iv = grid_iv[(grid_iv >= 0.0) & (grid_iv <= 4.4)]
            rows = []
            warm_u = np.zeros(params.sites)
            for rs in grid_iv:
                res = analyze_point(replace(params, r_s=float(rs)), eps,
                                    warm_u, rng, args.k_random)
                res["r_s"] = float(rs)
                warm_u = res.pop("warm_next")
                res.pop("u", None)
                rows.append(res)
                print(f"  [{name}] rs={rs:.3f} adm={res['admissible']} "
                      f"min={res['n_minima']} contact={res['n_ordering_contacts']} "
                      f"fail={res['n_failures']}", flush=True)
            # bisection of admissibility edges within this interval
            edges = []
            for i in range(len(rows) - 1):
                if rows[i]["admissible"] != rows[i + 1]["admissible"]:
                    edges.append((rows[i]["r_s"], rows[i + 1]["r_s"]))
            bounds = {}
            for e_a, e_b in edges:
                key = f"{e_a:.3f}/{e_b:.3f}"
                bounds[key] = bisect_edge(params, eps, e_a, e_b,
                                          np.zeros(params.sites), rng,
                                          args.k_random)
                print(f"    bisect {key} -> {bounds[key]:.4f}", flush=True)
            iv_out[name] = {"grid": rows, "bisected_edges": bounds,
                            "n_admissible": sum(r["admissible"] for r in rows)}
        per_eps[f"eps_frac_{eps_frac}"] = {"eps": eps, "intervals": iv_out}

    n_adm = {ef: sum(iv["n_admissible"]
                     for iv in d["intervals"].values())
             for ef, d in per_eps.items()}
    result = {
        "kind": "SSZ-FK ordering-constraint probe V2 (ChatGPT round-II task)",
        "scope": "EFFECTIVE_MODEL_ONLY; geometric ordering constraint "
                 "r_{i+1}-r_i>=eps replaces the per-site box; no "
                 "astrophysical or QNM claims",
        "motivation": "Box |u|<=0.45*dr is sufficient but not necessary for "
                      "ordering u_{i+1}-u_i > -dr; KKT boundary hit alone "
                      "does not prove a physical fold (counterexample "
                      "U=(u-mu)^2/2, |u|<=a).",
        "config": {"sites": args.sites, "eps_fracs": args.eps_fracs,
                   "k_random": args.k_random, "seed": args.seed,
                   "v0_bound": V0_BOUND, "gap_upper": GAP_UPPER,
                   "grad_tol": GRAD_TOL,
                   "intervals": intervals},
        "box_reference_summary": box_summary,
        "unconstrained_continuation": {
            "rows": unc_rows,
            "lambda_zero_crossings": crossings,
            "n_converged": sum(1 for r in unc_rows if r["converged"])},
        "per_eps": per_eps,
        "summary": {"admissible_points_per_eps": n_adm,
                    "unconstrained_lambda_crossings": len(crossings)},
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(f"\nwrote {args.output}")
    print(json.dumps(result["summary"], indent=2))


if __name__ == "__main__":
    main()
