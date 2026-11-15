# -*- coding: utf-8 -*-
"""Model loading and pathway-profile prediction for MASHres-NP."""
import json
import pickle
from importlib import resources

import numpy as np
import pandas as pd

from .featurize import featurize


class MASHresNP:
    """Structure-to-pathway model coupled with the in-house MASH disease vector.

    Usage
    -----
    >>> m = MASHresNP()
    >>> profile = m.predict_pathways("CC(=O)OC1=CC=CC=C1C(=O)O")   # pd.Series, 372 KEGG pathways
    >>> score = m.reversal_score("CC(=O)OC1=CC=CC=C1C(=O)O")       # float; higher = more MASH-reversing
    """

    def __init__(self):
        asset = resources.files("mashres_np") / "assets"
        pk = pickle.load(open(asset / "np_ridge_cellz_q50.pkl", "rb"))
        self.model, self.sx, self.sy = pk["model"], pk["x_scaler"], pk["y_scaler"]
        self.pathways = [p.replace("hsa", "map") for p in pk["pathways"]]
        self.pathway_names = json.load(open(asset / "kegg_pathway_names.json"))
        dvec = pd.read_csv(asset / "nash_disease_vector.csv", index_col=0).iloc[:, 0]
        dvec.index = [c.replace("hsa", "map") for c in dvec.index]
        self.disease_vector = dvec
        self.panel = self._build_panel()

    def _build_panel(self):
        """MASH-relevant pathway panel (keyword match on KEGG pathway names)."""
        kw = ["fatty acid", "lipid", "cholesterol", "steroid", "bile acid", "adipocyte",
              "insulin", "ppar", "ampk", "mtor", "diabetic", "tnf", "il-17", "nf-kappa",
              "nod-like", "toll-like", "chemokine", "cytokine", "inflammation", "complement",
              "ecm", "tgf", "glutathione", "oxidative phosphorylation", "glycolysis", "tca",
              "pentose", "gluconeogen", "arachidonic", "linoleic", "ether lipid",
              "glycerolipid", "glycerophospholipid", "sphingolipid", "peroxisome",
              "adipogenesis", "age-rage"]
        return [p for p in self.pathways
                if any(k in str(self.pathway_names.get(p, "")).lower() for k in kw)]

    def predict_pathways(self, smiles: str) -> pd.Series | None:
        """Predict the 372-pathway response profile for one compound."""
        v = featurize(smiles)
        if v is None:
            return None
        Xs = np.nan_to_num(self.sx.transform(v.reshape(1, -1)))
        pred = (self.model.predict(Xs) * self.sy.scale_ + self.sy.mean_)[0]
        return pd.Series(pred, index=self.pathways)

    def reversal_score(self, smiles: str) -> float | None:
        """MASH reversal score: negative Pearson correlation between the predicted
        pathway profile and the MASH disease vector on the shared panel."""
        prof = self.predict_pathways(smiles)
        if prof is None:
            return None
        common = [p for p in self.panel
                  if p in self.disease_vector.index and pd.notna(self.disease_vector[p])]
        dv = self.disease_vector[common].astype(float).values
        return float(-np.corrcoef(prof[common].values, dv)[0, 1])
