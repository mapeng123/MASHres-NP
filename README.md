# MASHres-NP

**MASH reversal screening of natural products** — a structure-to-pathway machine-learning pipeline that ranks natural compounds by their predicted ability to reverse the transcriptional pathway program of metabolic dysfunction-associated steatohepatitis (MASH).

## Overview

MASHres-NP couples two assets:

1. **Structure → pathway model.** A ridge-regression model trained on natural-product perturbation signatures from LINCS L1000 (cell-line z-score profiles, q50 aggregated), mapping a compound's Morgan fingerprints (2048 bits, r = 2) + MACCS keys (167) + full RDKit descriptor set (2432 features in total) to a **372-dimensional KEGG pathway response profile**.
2. **MASH disease vector.** An in-house KEGG pathway vector derived from differential analysis of human liver biopsies (GSE135251; NASH vs. healthy), capturing the transcriptional pathway program of the disease.

The **MASH reversal score** of a compound is the negative Pearson correlation between its predicted pathway profile and the MASH disease vector, computed on a curated panel of 49 MASH-relevant pathways (lipid metabolism, insulin/AMPK/PPAR signaling, inflammation, ECM/fibrosis, oxidative metabolism). A higher score indicates a stronger predicted tendency to reverse the MASH transcriptional program.

The pipeline rediscovered **oxypaeoniflorin (OPA)** — a monoterpene glycoside from *Paeonia lactiflora* — as a top candidate, a finding validated in PA-induced AML12 hepatocytes and in a high-fat-diet mouse model of MASH (see the accompanying manuscript).

## Installation

```bash
git clone https://github.com/mapeng123/MASHres-NP.git
cd MASHres-NP
pip install -e .
```

Requires Python ≥ 3.10. Dependencies: `rdkit`, `scikit-learn` (≥ 1.9), `pandas`, `numpy`.

## Quick start

### 1. Screen your own compound library

Prepare a CSV with at least a `smiles` column (a `name` column is optional):

```csv
name,smiles
my_compound,CC(=O)OC1=CC=CC=C1C(=O)O
```

Then run:

```bash
mashres-np screen --input my_library.csv --output my_library_ranked.csv
```

The output is ranked by `mash_reversal_score` (descending); unparseable SMILES are dropped.

### 2. Python API

```python
from mashres_np import MASHresNP, screen_library
import pandas as pd

m = MASHresNP()

# single-compound reversal score
score = m.reversal_score("CC(=O)OC1=CC=CC=C1C(=O)O")

# full 372-pathway predicted profile
profile = m.predict_pathways("CC(=O)OC1=CC=CC=C1C(=O)O")

# batch screening
lib = pd.DataFrame({"name": ["aspirin"], "smiles": ["CC(=O)OC1=CC=CC=C1C(=O)O"]})
ranked = screen_library(lib, model=m)
```

### 3. Reproduce the demo

```bash
mashres-np screen --input examples/demo_compounds.csv --output examples/demo_ranked.csv
```

Expected output (rank / score):

| rank | compound | score |
|------|----------|-------|
| 1 | baicalin | ≈ 0.494 |
| 2 | aspirin | ≈ 0.448 |
| 3 | **oxypaeoniflorin (OPA)** | **≈ 0.400** |
| 4 | albiflorin | ≈ 0.400 |
| 5 | paeoniflorin | ≈ 0.388 |

OPA reproduces a reversal score of ≈ 0.3999, consistent with the value reported in the manuscript.

## Repository layout

```
MASHres-NP/
├── mashres_np/
│   ├── __init__.py
│   ├── featurize.py        # Morgan FP + MACCS + RDKit descriptors (2432-d)
│   ├── model.py            # MASHresNP: model loading, pathway prediction, scoring
│   ├── screen.py           # batch library screening
│   ├── cli.py              # `mashres-np` command-line interface
│   └── assets/
│       ├── np_ridge_cellz_q50.pkl      # trained ridge model + scalers
│       ├── nash_disease_vector.csv     # MASH disease vector (KEGG pathways)
│       └── kegg_pathway_names.json     # pathway ID -> name map
├── examples/
│   └── demo_compounds.csv
├── tests/
│   └── test_reproduce.py
├── pyproject.toml
├── requirements.txt
└── LICENSE
```

## Method notes

- **Feature pipeline** (must not be modified): Morgan fingerprint (radius 2, 2048 bits) + MACCS keys (167 bits) + the full sorted RDKit descriptor list, concatenated in this order; features are standardized with the training-fitted scaler (NaNs are zero-filled after scaling).
- **Prediction**: ridge regression, predictions de-standardized with the pathway-scaler fitted on the LINCS training set.
- **Scoring panel**: 49 KEGG pathways whose names match MASH-relevant keywords (fatty acid/lipid metabolism, insulin/PPAR/AMPK, inflammatory signaling, ECM/TGF-β, oxidative phosphorylation, arachidonic/linoleic acid metabolism, …), intersected with non-missing entries of the disease vector.
- **Reversal score**: `-Pearson r` between the predicted pathway profile and the MASH disease vector on the shared panel. Scores are relative, not absolute; use them to rank compounds within one library.

## Citation

If you use MASHres-NP, please cite:

> [Manuscript under review] MASHres-NP: Machine Learning-Guided Screening of Natural Products Identifies Oxypaeoniflorin as a MASH-Reversing Agent.

## License

MIT (see `LICENSE`). The bundled model and disease vector are released for academic use.
