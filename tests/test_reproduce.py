# -*- coding: utf-8 -*-
"""Reproduction test: the bundled model must reproduce the reported OPA score
and the expected demo ranking.

Run either with pytest (`pytest tests/`) or standalone from the repo root:
    python tests/test_reproduce.py
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd

from mashres_np import MASHresNP, screen_library

OPA = ("C[C@@]12C[C@@]3(O)O[C@@H](O1)[C@]1(COC(=O)C4=CC=C(O)C=C4)"
       "[C@H]3CC12O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O")


def test_opa_score_reproduced():
    m = MASHresNP()
    s = m.reversal_score(OPA)
    assert s is not None and math.isclose(s, 0.3999, abs_tol=1e-3), f"OPA score {s}"


def test_panel_size():
    m = MASHresNP()
    assert len(m.pathways) == 372
    assert len(m.panel) == 49


def test_demo_ranking():
    m = MASHresNP()
    lib = pd.read_csv(os.path.join(os.path.dirname(__file__), "..", "examples", "demo_compounds.csv"))
    out = screen_library(lib, model=m)
    assert out.iloc[0]["name"] == "baicalin"
    assert out.iloc[-1]["name"] == "paeoniflorin"
    assert "oxypaeoniflorin" in set(out["name"])


if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")
    for fn in (test_opa_score_reproduced, test_panel_size, test_demo_ranking):
        fn()
        print(f"PASS  {fn.__name__}")
    print("All reproduction tests passed.")
