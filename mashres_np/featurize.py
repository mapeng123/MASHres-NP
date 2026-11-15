# -*- coding: utf-8 -*-
"""Molecular featurization for MASHres-NP.

Exactly replicates the feature pipeline used to train the bundled model:
Morgan FP (radius 2, 2048 bits) + MACCS keys (167) + full RDKit descriptor set
(sorted descList, float64). Feature order must NOT be changed.
"""
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, MACCSkeys, Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

RDLogger.DisableLog("rdApp.*")

_DESC_NAMES = sorted(n for n, _ in Descriptors.descList)
_CALC = MoleculeDescriptors.MolecularDescriptorCalculator(_DESC_NAMES)
N_FEATURES = 2048 + 167 + len(_DESC_NAMES)


def featurize(smiles: str) -> np.ndarray | None:
    """Return the 2432-d feature vector for one SMILES, or None if unparseable."""
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return None
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
    arr = np.zeros(2048, dtype=np.float32)
    arr[list(fp.GetOnBits())] = 1.0
    mac = MACCSkeys.GenMACCSKeys(mol)
    marr = np.zeros(167, dtype=np.float32)
    marr[list(mac.GetOnBits())] = 1.0
    dv = np.asarray(_CALC.CalcDescriptors(mol), dtype=np.float64)
    dv = np.nan_to_num(dv, nan=np.nan, posinf=np.nan, neginf=np.nan)
    return np.concatenate([arr, marr, dv])
