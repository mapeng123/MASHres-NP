# -*- coding: utf-8 -*-
"""Batch screening: score a SMILES library and return a ranked table."""
from __future__ import annotations

import pandas as pd

from .model import MASHresNP


def screen_library(df: pd.DataFrame,
                   smiles_col: str = "smiles",
                   name_col: str | None = "name",
                   model: MASHresNP | None = None) -> pd.DataFrame:
    """Score every compound in ``df`` and return a table ranked by MASH
    reversal score (descending).

    Parameters
    ----------
    df : DataFrame containing at least a SMILES column.
    smiles_col : name of the SMILES column (default "smiles").
    name_col : name of the compound-name column, or None.
    model : optional pre-loaded MASHresNP instance (loaded once if omitted).

    Returns
    -------
    DataFrame with columns [rank, name?, smiles, mash_reversal_score],
    unparseable SMILES dropped with score NaN removed.
    """
    if smiles_col not in df.columns:
        raise ValueError(f"column '{smiles_col}' not found in input table")
    m = model or MASHresNP()

    records = []
    for _, row in df.iterrows():
        smi = str(row[smiles_col]).strip()
        if not smi or smi.lower() == "nan":
            continue
        score = m.reversal_score(smi)
        rec = {"smiles": smi, "mash_reversal_score": score}
        if name_col and name_col in df.columns:
            rec["name"] = row[name_col]
        records.append(rec)

    out = pd.DataFrame(records).dropna(subset=["mash_reversal_score"])
    out = out.sort_values("mash_reversal_score", ascending=False).reset_index(drop=True)
    out.insert(0, "rank", range(1, len(out) + 1))
    cols = ["rank"] + (["name"] if name_col and "name" in out.columns else []) + \
           ["smiles", "mash_reversal_score"]
    return out[cols]
