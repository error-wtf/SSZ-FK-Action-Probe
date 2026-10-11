#!/usr/bin/env python3
"""Demo: byte-level reproducibility of the whole package pipeline.

Re-runs the continuation step into a temp dir and compares sha256 with the
shipped artifact SSZ_FK_CURVATURE_CONTINUATION_V1.json. A mismatch is a
FAIL (determinism broken). The full pipeline (3 pilots + 12 tests + 89-step
continuation + saddle zones) was verified byteidentical on 2026-10-11:
  pilot V1  8d88161a…   N60 fac13449…
  continuation de114e1d…   saddle 33a3e921…

Run:  /home/error/physics/venv/bin/python demos/demo_reproducibility.py
Time: ~2-3 min (89 continuation steps).
"""
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

PKG = Path(__file__).resolve().parent.parent
PY = sys.executable

def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    print(__doc__)
    ref = PKG / "SSZ_FK_CURVATURE_CONTINUATION_V1.json"
    out = Path(tempfile.mkdtemp(prefix="demo_repro_")) / "cont.json"
    print("re-running continuation (89 steps) ...")
    r = subprocess.run([PY, str(PKG / "tools/curvature_continuation_v1.py"),
                        "--steps", "89", "--output", str(out)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print("CONTINUATION FAILED:", r.stderr[-800:]); sys.exit(1)
    h_ref, h_new = sha(ref), sha(out)
    print(f"shipped : {h_ref}")
    print(f"re-run  : {h_new}")
    ok = h_ref == h_new
    print("REPRODUCIBLE: YES" if ok else "REPRODUCIBLE: NO  <-- INVESTIGATE")
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
