# -*- coding: utf-8 -*-
"""MASHres-NP: MASH reversal screening of natural products.

A structure-to-pathway ridge model (trained on LINCS L1000 natural-product
perturbation profiles) coupled with an in-house MASH disease vector derived
from human liver biopsies (GSE135251), ranking natural compounds by their
predicted MASH-reversal potential.
"""
from .model import MASHresNP
from .screen import screen_library

__version__ = "0.1.0"
__all__ = ["MASHresNP", "screen_library", "__version__"]
