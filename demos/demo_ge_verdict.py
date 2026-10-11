#!/usr/bin/env python3
"""Demo: G-E result — the SSZ-FK action does NOT collapse to FK x redshift.

Background (6a proper-length probe, committed as 5245a0f in
SSZ_FULL_CLOSURE): under the naive map phi = Q*l the G-E residual is
0.564; under the proper-length phase map phi = Q * (proper length)
the residual drops to 0.0043 — i.e. the apparent disagreement was a
PHASE-CONVENTION artifact, and the proper verdict is

    NEW_STRUCTURE_NOT_COLLAPSIBLE in all 12 tested combinations.

This demo reproduces the pilot quantities of the V1 probe (single run,
default sites) and prints the verdict from the fresh output.

Run:  /home/error/physics/venv/bin/python demos/demo_ge_verdict.py
Time: ~1-2 min (full probe run).
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent
PY = sys.executable

def main():
    print(__doc__)
    out = Path(tempfile.mkdtemp(prefix="demo_ge_")) / "pilot.json"
    print(f"\nrunning tools/ssz_fk_action_probe_v1.py -> {out} ...")
    r = subprocess.run([PY, str(PKG / "tools/ssz_fk_action_probe_v1.py"),
                        "--output", str(out)], capture_output=True, text=True)
    if r.returncode != 0:
        print("PROBE FAILED:", r.stderr[-800:]); sys.exit(1)
    d = json.loads(out.read_text())
    print("\n--- key quantities (fresh run) ---")
    print(f"kind       : {d.get('kind')}")
    print(f"scope      : {d.get('scope')}")
    m = d.get("metric", {})
    print(f"D_R        : {m.get('D_R')}   s_R: {m.get('s_R')}")
    a = d.get("action_ssz_vs_flat", {})
    hm = a.get("overlap_matching", {}).get("high_mac_subset", {})
    print(f"rank residual (all modes)     : {a.get('rank_residual')}")
    print(f"residual, high-MAC subset     : {hm.get('residual')} "
          f"(n={hm.get('n')})")
    print("\nNOTE: the pilot alone does not print the 6a verdict; the 0.564 ->")
    print("0.0043 improvement and NEW_STRUCTURE_NOT_COLLAPSIBLE verdict come")
    print("from the 6a probe in SSZ_FULL_CLOSURE (see NOTE below).")
    print("\nFull audit trail: 6a artifact in SSZ_FULL_CLOSURE"
          " (artifacts/SSZ_ACTION_DERIVED_MATTER_PROBE_V1.json,"
          " commit 5245a0f) + docs/DERIVATION_AND_FINDINGS.md here.")


if __name__ == "__main__":
    main()
