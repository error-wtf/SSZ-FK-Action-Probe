#!/usr/bin/env python3
"""Figures for the SSZ-FK proper-length research package (V1).

Generates static PNGs + animated GIFs from the REAL artifacts of this
package (no synthetic data):

  figures/fig1_regime_map.png        static  — energy + stability vs r_s,
                                               admissibility from multistart,
                                               bisected fold borders
  figures/fig2_mode_tracking.png     static  — MAC/subspace matching
                                               diagnostics vs r_s
  figures/anim_chain_equilibrium.gif animated — ground-state chain
                                               configuration along the
                                               r_s continuation (re-solved
                                               live, warm-started)
  figures/anim_admissibility_scan.gif animated — the admissibility map
                                               being built up during the
                                               scan (coarse map -> gap
                                               sweeps -> bisected borders)

Usage:
  ../venv python tools/make_figures_v1.py [--outdir figures] [--fps 6]

Scope: EFFECTIVE_MODEL_ONLY — no astrophysical or QNM claims.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import replace
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

PKG = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PKG / "tools"))

from ssz_fk_action_probe_v1 import Parameters, SSZFKAction  # noqa: E402
from curvature_continuation_v1 import warm_equilibrium  # noqa: E402


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def load_artifacts():
    cc = json.loads((PKG / "SSZ_FK_CURVATURE_CONTINUATION_V1.json").read_text())
    sz = json.loads((PKG / "SSZ_FK_SADDLE_ZONES_V1.json").read_text())
    return cc, sz


def merged_admissibility(sz: dict):
    """Coarse map + gap sweeps -> (rs, admissible) sorted by rs."""
    pts = []
    for rec in sz["coarse_map"]:
        pts.append((rec["r_s"], rec["admissible"]))
    for sweep in sz["gap_sweeps"].values():
        for rec in sweep:
            pts.append((rec["r_s"], rec["admissible"]))
    pts.sort(key=lambda t: t[0])
    return pts


def fig1_regime_map(cc: dict, sz: dict, out: Path):
    rs = np.array([b["r_s"] for b in cc["branch"]])
    energy = np.array([b["energy"] for b in cc["branch"]])
    mineig = np.array([b["min_hessian_eig"] for b in cc["branch"]])
    adm = merged_admissibility(sz)
    rs_a = np.array([p[0] for p in adm])
    ok_a = np.array([p[1] for p in adm])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7.5), sharex=True,
                                   gridspec_kw={"height_ratios": [3, 2]})

    # -- top: energy + admissibility band
    ax1.scatter(rs_a[ok_a], np.full(ok_a.sum(), 0.92), s=14, color="tab:green",
                marker="s", label="admissible (multistart, 9 starts)")
    ax1.scatter(rs_a[~ok_a], np.full((~ok_a).sum(), 0.92), s=14, color="tab:red",
                marker="s", label="CONSTRAINT_BOUNDARY (no interior GG)")
    ax1.plot(rs, energy, "-", color="tab:blue", lw=1.4,
             label="ground-state energy (continuation)")
    ax1.set_ylabel("E(u*)  [a.u.]")
    ax1.set_title("SSZ-FK ground state along the curvature continuation\n"
                  "gaps = CONSTRAINT_BOUNDARY (fold at the ordering domain); "
                  "no saddles, no hysteresis")
    # bisected borders
    for key, b in sz["boundaries_bisected"].items():
        ax1.axvline(b, color="k", ls=":", lw=0.9)
        ax1.annotate(f"{b:.4f}", xy=(b, ax1.get_ylim()[1]*0.02),
                     fontsize=6.5, rotation=90, va="bottom", ha="right")
    ax1.legend(loc="upper left", fontsize=8)
    ax1.grid(alpha=0.3)

    # -- bottom: Hessian min eig (stability)
    ax2.axhline(0.0, color="k", lw=0.8)
    ax2.plot(rs, mineig, "-", color="tab:orange", lw=1.4,
             label="min Hessian eigenvalue")
    ax2.scatter(rs_a[~ok_a], np.full((~ok_a).sum(), 0.0), s=14,
                color="tab:red", marker="s")
    ax2.set_xlabel("r_s  (curvature parameter, background family scan)")
    ax2.set_ylabel("lambda_min(H)")
    ax2.set_yscale("symlog", linthresh=1e-3)
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fig2_mode_tracking(cc: dict, out: Path):
    rs = np.array([b["r_s"] for b in cc["branch"]])
    alpha_mac = np.array([b.get("mac", {}).get("alpha", np.nan)
                          for b in cc["branch"]])
    mac_min = np.array([b.get("mac", {}).get("mac_min", np.nan)
                        for b in cc["branch"]])
    resid_sub = np.array([b.get("subspace", {}).get("max_rel_resid", np.nan)
                          for b in cc["branch"]])
    n_match = int(np.isfinite(alpha_mac).sum())
    print(f"[fig2] matched steps with MAC diagnostics: {n_match}/{len(rs)}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4), sharex=True)
    ax1.plot(rs, alpha_mac, "-", lw=1.3, color="tab:purple",
             label="MAC alpha (rank-matched)")
    ax1.plot(rs, mac_min, "-", lw=1.0, color="tab:cyan",
             label="MAC min (top-20)")
    ax1.set_ylim(-0.05, 1.05)
    ax1.set_xlabel("r_s"); ax1.set_ylabel("MAC")
    ax1.set_title("mode-match diagnostics")
    ax1.legend(fontsize=8); ax1.grid(alpha=0.3)

    ax2.plot(rs, resid_sub, "-", lw=1.3, color="tab:brown",
             label="subspace max rel. residual")
    ax2.set_yscale("log")
    ax2.set_xlabel("r_s"); ax2.set_ylabel("max rel. resid")
    ax2.set_title("subspace residual (eigenpair quality)")
    ax2.legend(fontsize=8); ax2.grid(alpha=0.3)
    fig.suptitle("Mode identity across regimes: jumps in MAC/residual at the "
                 "fold borders\n(rank- and MAC-matching see the same jumps -> "
                 "different mode family after each gap)", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(out, dpi=150)
    plt.close(fig)


def anim_chain_equilibrium(out: Path, fps: int = 6):
    """Re-solve the warm continuation live and animate the chain shape."""
    p0 = Parameters(sites=40)
    steps = np.linspace(0.0, 4.4, 89)

    # pre-solve all equilibria (warm continuation, identical to gate 1)
    frames = []
    u_prev = np.zeros(40)
    for rs in steps:
        p = replace(p0, r_s=float(rs))
        mdl = SSZFKAction(p, flat=False)
        eq = warm_equilibrium(mdl, u_prev) if rs > 0 else \
            warm_equilibrium(mdl, np.zeros(40))
        u_prev = eq["u"]
        r = mdl.position(eq["u"])
        frames.append((float(rs), r, float(eq["energy"]),
                       float(eq["min_hessian_eig"]),
                       bool(eq["stationary"])))

    fig, ax = plt.subplots(figsize=(8, 5))
    line, = ax.plot([], [], "o-", ms=3, lw=1.0, color="tab:blue")
    txt = ax.text(0.02, 0.95, "", transform=ax.transAxes, fontsize=10,
                  va="top", family="monospace")

    rs_all = [f[0] for f in frames]
    r_all = np.concatenate([f[1] for f in frames])
    ax.set_xlim(min(r_all) - 0.2, max(r_all) + 0.2)
    ax.set_ylim(-0.6, 0.6)
    ax.set_xlabel("site index k (radial position in legend: proper-length map)")
    ax.set_ylabel("chain displacement (exaggerated)")
    ax.grid(alpha=0.3)

    def update(i):
        rs, r, e, me, stat = frames[i]
        # deviation of each site from a homogeneous spacing reference:
        #   d_k = r_k - r_0 - k * mean_spacing
        dev = r - r[0] - np.arange(len(r)) * np.mean(np.diff(r))
        line.set_data(np.arange(len(r)), dev)
        txt.set_text(f"r_s = {rs:.2f}\nE = {e:.5f}\n"
                     f"lambda_min(H) = {me:.4f}\nstationary = {stat}")
        return line, txt

    ax.set_ylabel("deviation from homogeneous spacing")
    anim = FuncAnimation(fig, update, frames=len(frames), interval=1000 // fps)
    anim.save(out, writer=PillowWriter(fps=fps))
    plt.close(fig)
    return frames


def anim_admissibility_scan(sz: dict, out: Path, fps: int = 6):
    pts = merged_admissibility(sz)
    borders = sorted(sz["boundaries_bisected"].values())

    fig, ax = plt.subplots(figsize=(9, 2.6))
    ok_pts = [(r, a) for r, a in pts]
    sc_ok = ax.scatter([], [], s=26, marker="s", color="tab:green",
                       label="admissible")
    sc_bad = ax.scatter([], [], s=26, marker="s", color="tab:red",
                        label="CONSTRAINT_BOUNDARY")
    txt = ax.text(0.01, 1.25, "", transform=ax.transAxes, fontsize=9,
                  family="monospace")
    ax.set_xlim(-0.08, 4.55)
    ax.set_ylim(-0.6, 1.6)
    ax.set_yticks([])
    ax.set_xlabel("r_s")
    ax.legend(loc="upper center", ncol=2, fontsize=8, frameon=False)
    ax.set_title("admissibility map build-up (coarse 0.1 -> gap 0.025 -> "
                 "bisection ±0.005)", fontsize=9)

    n_coarse = len(sz["coarse_map"])
    g_ok_x, g_ok_y, g_bad_x, g_bad_y = [], [], [], []

    def update(i):
        r, a = ok_pts[i]
        (g_ok_x if a else g_bad_x).append(r)
        (g_ok_y if a else g_bad_y).append(0.0)
        phase = ("coarse map" if i < n_coarse else
                 "gap sweeps (dr=0.025)" if i < len(pts) - 0 else "?")
        if a:
            sc_ok.set_offsets(np.c_[g_ok_x, g_ok_y])
        else:
            sc_bad.set_offsets(np.c_[g_bad_x, g_bad_y])
        txt.set_text(f"[{i+1}/{len(pts)}] r_s={r:.3f} -> "
                     f"{'MINIMUM' if a else 'no interior GG'}  |  {phase}")
        return sc_ok, sc_bad, txt

    anim = FuncAnimation(fig, update, frames=len(pts), interval=1000 // fps)
    anim.save(out, writer=PillowWriter(fps=fps))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", type=Path, default=PKG / "figures")
    ap.add_argument("--fps", type=int, default=6)
    a = ap.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)

    cc, sz = load_artifacts()
    made = []

    f1 = a.outdir / "fig1_regime_map.png"
    fig1_regime_map(cc, sz, f1); made.append(f1)
    f2 = a.outdir / "fig2_mode_tracking.png"
    fig2_mode_tracking(cc, f2); made.append(f2)
    f3 = a.outdir / "anim_chain_equilibrium.gif"
    anim_chain_equilibrium(f3, fps=a.fps); made.append(f3)
    f4 = a.outdir / "anim_admissibility_scan.gif"
    anim_admissibility_scan(sz, f4, fps=a.fps); made.append(f4)

    print("== figures written ==")
    for p in made:
        print(f"{p.name}: {p.stat().st_size} bytes  sha256={_sha256(p)}")


if __name__ == "__main__":
    main()
