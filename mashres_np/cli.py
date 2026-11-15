# -*- coding: utf-8 -*-
"""Command-line interface for MASHres-NP."""
from __future__ import annotations

import argparse
import sys

import pandas as pd

from .model import MASHresNP
from .screen import screen_library


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="mashres-np",
        description="MASHres-NP: rank natural products by predicted MASH-reversal potential.")
    sub = ap.add_subparsers(dest="cmd", required=True)

    p_scr = sub.add_parser("screen", help="Screen a compound library (CSV with a SMILES column).")
    p_scr.add_argument("--input", required=True, help="Input CSV file.")
    p_scr.add_argument("--output", required=True, help="Output CSV file (ranked).")
    p_scr.add_argument("--smiles-col", default="smiles", help="SMILES column name (default: smiles).")
    p_scr.add_argument("--name-col", default="name", help="Name column name (default: name; use '' to disable).")

    p_pred = sub.add_parser("predict", help="Predict the 372-pathway profile of one compound.")
    p_pred.add_argument("--smiles", required=True, help="SMILES string.")
    p_pred.add_argument("--output", default=None, help="Optional output CSV (pathway, predicted_response).")

    args = ap.parse_args(argv)

    if args.cmd == "screen":
        df = pd.read_csv(args.input)
        name_col = args.name_col or None
        out = screen_library(df, smiles_col=args.smiles_col, name_col=name_col)
        out.to_csv(args.output, index=False)
        print(f"[MASHres-NP] screened {len(df)} compounds -> {len(out)} scored, "
              f"ranked table written to {args.output}", file=sys.stderr)
        print(out.head(10).to_string(index=False))
    elif args.cmd == "predict":
        m = MASHresNP()
        prof = m.predict_pathways(args.smiles)
        if prof is None:
            sys.exit("[MASHres-NP] could not parse SMILES.")
        df = pd.DataFrame({"pathway": prof.index, "predicted_response": prof.values})
        df["pathway_name"] = df["pathway"].map(m.pathway_names)
        if args.output:
            df.to_csv(args.output, index=False)
            print(f"[MASHres-NP] pathway profile written to {args.output}", file=sys.stderr)
        else:
            print(df.to_string(index=False))


if __name__ == "__main__":
    main()
