#!/usr/bin/env python3
"""Demo: gap-width systematics — fold borders narrow with curvature.

The three inadmissibility gaps have bisected borders

    gap 1:  0.2795 -> 0.8998   (width 0.6203)
    gap 2:  1.8580 -> 2.3498   (width 0.4918)
    gap 3:  2.6057 -> 2.9502   (width 0.3445)

i.e. the gap widths DECREASE with increasing r_s (curvature). This demo
prints the numbers from the real artifact and fits width vs centre to make
the trend explicit. It also states the open physics question: is this a
property of the model (physics) or of the box size 0.45*dr (numerics)?
The box-sensitivity study (0.30/0.45/0.49*dr) is the designated test.

Run:  /home/error/physics/venv/bin/python demos/demo_gap_systematics.py
Time: <1 s (reads the artifact only).
"""
import json
from pathlib import Path

import numpy as np

PKG = Path(__file__).resolve().parent.parent

def main():
    sz = json.loads((PKG / "SSZ_FK_SADDLE_ZONES_V1.json").read_text())
    print(__doc__)
    print("--- bisected borders (from SSZ_FK_SADDLE_ZONES_V1.json) ---")
    items = sorted(sz["boundaries_bisected"].items())
    pairs = [(items[i][1], items[i + 1][1]) for i in range(0, len(items), 2)]
    print(f"{'gap':<6}{'lo':>10}{'hi':>10}{'width':>10}{'centre':>10}")
    widths, centres = [], []
    for i, (lo, hi) in enumerate(pairs, 1):
        w = hi - lo
        widths.append(w); centres.append((lo + hi) / 2)
        print(f"{i:<6}{lo:>10.4f}{hi:>10.4f}{w:>10.4f}{centres[-1]:>10.4f}")

    # linear trend width(centre)
    A = np.vstack([centres, np.ones(3)]).T
    slope, intercept = np.linalg.lstsq(A, np.array(widths), rcond=None)[0]
    resid = np.array(widths) - (slope * np.array(centres) + intercept)
    print(f"\nlinear fit: width ~ {slope:.4f} * centre + {intercept:.4f}")
    print(f"fit residuals: {np.array2string(resid, precision=5)}")
    print("\nInterpretation: widths shrink monotonically with curvature")
    print("(0.62 -> 0.49 -> 0.34). OPEN: physics vs box artefact — decide")
    print("via box sensitivity 0.30/0.45/0.49*dr (priority 1, not yet run).")
    print("IMPORTANT scope note: the scan is a parametrised background")
    print("family, NOT a time evolution (ChatGPT round-I correction).")


if __name__ == "__main__":
    main()
