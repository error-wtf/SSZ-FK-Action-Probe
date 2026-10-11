#!/usr/bin/env python3
"""Figures for the ordering-constraint probe (real artifact data only).

fig3_ordering_constraint.png
  (a) unconstrained lambda_min over r_s (0..4.4) with convergence markers
      and the lambda=0 line;
  (b) admissibility strips for the three former box gaps under the true
      ordering constraint (eps = 0.05/0.2/0.5 dr) vs. the box-constrained
      classification (red) from SSZ_FK_SADDLE_ZONES_V1.json;
  (c) minimal radial gap of the unconstrained solutions (ordering check).

anim_unconstrained_continuation.gif
  Chain configuration u_i(r_s) along the unconstrained Newton path
  (recomputed live, same discipline as anim_chain_equilibrium.gif);
  frames where Newton stalls are labelled explicitly as non-stationary.
"""
from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

import sys
sys.path.insert(0, str(Path(__file__).parent))
from ordering_constraint_probe_v1 import unconstrained_newton
from ssz_fk_action_probe_v1 import Parameters, SSZFKAction

ROOT = Path(__file__).resolve().parent.parent
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)


def main():
    art = json.loads((ROOT / "SSZ_FK_ORDERING_CONSTRAINT_V1.json").read_text())
    box = json.loads((ROOT / "SSZ_FK_SADDLE_ZONES_V1.json").read_text())

    unc = art["unconstrained_continuation"]["rows"]
    rs_u = np.array([r["r_s"] for r in unc])
    lam = np.array([r["min_hessian_eig"] for r in unc])
    conv = np.array([r["converged"] for r in unc])
    mingap = np.array([r["min_gap"] for r in unc])

    eps_keys = sorted(art["per_eps"].keys(),
                      key=lambda k: float(k.split("_")[-1]))
    gaps = {"gap1": (0.2795, 0.8998), "gap2": (1.8580, 2.3498),
            "gap3": (2.6057, 2.9502)}

    # box classification inside the gaps (gap_sweeps: dict gap -> list)
    box_cls = {}
    for gname, rows in box.get("gap_sweeps", {}).items():
        for row in rows:
            r = row.get("r_s", row.get("rs"))
            if r is None:
                continue
            box_cls.setdefault(gname, {})[r] = row.get(
                "classification", row.get("kind"))

    fig, axes = plt.subplots(3, 1, figsize=(11, 12))

    ax = axes[0]
    ax.plot(rs_u[conv], lam[conv], "o-", ms=4, color="tab:blue",
            label=r"$\lambda_{\min}$ (Newton converged)")
    ax.plot(rs_u[~conv], lam[~conv], "x", ms=6, color="tab:red",
            label="Newton stall (no stationary point on this path)")
    ax.axhline(0.0, color="k", ls="--", lw=0.8)
    ax.set_ylabel(r"$\lambda_{\min}(H)$")
    ax.set_yscale("symlog", linthresh=1e-2)
    ax.set_title("Unconstrained Newton continuation (no ordering box)")
    ax.legend(loc="lower left", fontsize=9)

    ax = axes[1]
    yticks, ylabels = [], []
    for yi, gname in enumerate(gaps):
        a, b = gaps[gname]
        ax.add_patch(plt.Rectangle((a, yi - 0.33), b - a, 0.31,
                                   color="tab:red", alpha=0.35))
        for zi, ek in enumerate(eps_keys):
            iv = art["per_eps"][ek]["intervals"][gname]
            rs_i = [r["r_s"] for r in iv["grid"]]
            y = yi + 0.045 * (zi + 1)
            adm = [r["admissible"] for r in iv["grid"]]
            ax.scatter(rs_i, [y] * len(rs_i), s=22,
                       color=["tab:green" if x else "k" for x in adm],
                       marker="s", zorder=5)
        yticks.append(yi)
        ylabels.append(f"{gname} (box gap {a}–{b})")
    ax.set_yticks(yticks, ylabels)
    ax.set_ylim(-0.5, len(gaps) - 0.4)
    ax.set_xlabel(r"$r_s$")
    ax.set_title("Former box gaps under ordering constraint "
                 r"$r_{i+1}-r_i\geq\varepsilon$: red = box-inadmissible, "
                 "green = admissible "
                 r"($\varepsilon=0.05/0.2/0.5\,\Delta r$)")

    ax = axes[2]
    ax.plot(rs_u, mingap, "s-", ms=3, color="tab:purple")
    ax.set_xlabel(r"$r_s$")
    ax.set_ylabel(r"$\min_i (r_{i+1}-r_i)$")
    ax.set_title("Ordering check of unconstrained solutions "
                 f"(global min {mingap.min():.4f} > 0)")
    fig.suptitle("SSZ-FK ordering-constraint probe V1: the former 'gaps' are "
                 "box artefacts, not folds", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(FIG / "fig3_ordering_constraint.png", dpi=130)
    plt.close(fig)
    print("fig3_ordering_constraint.png written")

    # ---- animation: unconstrained chain over r_s (live recompute) ----
    p = Parameters(sites=40)
    q = np.arange(p.sites, dtype=float) * p.dr
    n_frames = 60
    rs_frames = np.linspace(0.0, 4.4, n_frames)
    warm = unconstrained_newton(SSZFKAction(replace(p, r_s=0.0), flat=False),
                                np.zeros(p.sites))

    fig2, ax2 = plt.subplots(figsize=(10, 4.2))
    rs_txt = ax2.text(0.02, 0.95, "", transform=ax2.transAxes, va="top",
                      fontsize=10)
    line, = ax2.plot([], [], lw=1.4)
    ax2.set_xlabel("site index k")
    ax2.set_ylabel(r"$u_k$ (deviation from homogeneous lattice)")
    ax2.set_xlim(0, p.sites - 1)

    def update(fi):
        nonlocal warm
        rs = float(rs_frames[fi])
        res = unconstrained_newton(
            SSZFKAction(replace(p, r_s=rs), flat=False),
            np.array(warm["u"]))
        u = np.array(res["u"])
        if res["converged"]:
            warm = res
            rs_txt.set_text(f"$r_s$={rs:.2f}  STATIONARY  "
                            f"$\\lambda_\\min$={res['min_hessian_eig']:.3e}  "
                            f"min gap={res['min_gap']:.4f}")
        else:
            rs_txt.set_text(f"$r_s$={rs:.2f}  NEWTON STALL — no stationary "
                            "point reached on this path (non-stationary "
                            "configuration shown)")
        line.set_data(np.arange(p.sites), u)
        ax2.set_ylim(-0.6 * p.dr, 0.6 * p.dr)
        return line, rs_txt

    meta = {"title": "unconstrained continuation of the SSZ-FK chain"}
    anim = FuncAnimation(fig2, update, frames=n_frames, interval=167)
    anim.save(FIG / "anim_unconstrained_continuation.gif",
              writer=PillowWriter(fps=6), dpi=100)
    plt.close(fig2)
    print("anim_unconstrained_continuation.gif written")


if __name__ == "__main__":
    main()
