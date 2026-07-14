"""Compute scores with fitdock
"""
"""Run an external scoring subprocess

Run external process: provide specific command line parameters when needed
pass on the SMILES as a series of strings at the end.
"""

#from __future__ import annotations

__all__ = ["Scorefitdock"]

import os
import logging
import copy

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np

from .component_results import ComponentResults
from .run_program import run_command
from .add_tag import add_tag
from .FitDock import fitscore
import sys


logger = logging.getLogger(__name__)


@add_tag("__parameters")
@dataclass
class Parameters:
    """Parameters for the scoring component

    Note that all parameters are always lists because components can have
    multiple endpoints and so all the parameters from each endpoint is
    collected into a list.  This is also true in cases where there is only one
    endpoint.
    """

    protein_pdb: List[str] 
    output_dir: List[str]
    ref_ligand_mol2: List[str]
    ddg: List[float]
    rmsd: List[float]


@add_tag("__component")
class Scorefitdock:
    """Run docking with fitdock and score by pbcnet2.0
    """

    def __init__(self, params: Parameters):

        self._internal_step = 0
        
        self.protein_pdb = params.protein_pdb[0]
        self.output_dir = params.output_dir[0]
        self.ref_ligand_mol2 = params.ref_ligand_mol2[0]
        self.ddg = params.ddg[0]
        self.rmsd = params.rmsd[0]

    def __call__(self, smilies: List[str]) -> np.array:
        scores = []
        #
        try:
            #set output dir
            if not os.path.exists(self.output_dir):
                os.makedirs(self.output_dir)

            #run fitdock to get ddG
            output_subdir = os.path.join(self.output_dir,f"epoch_{self._internal_step}")
            ddg_score = fitscore(output_dir=output_subdir,
                    smiles=smilies,
                    protein_pdb=self.protein_pdb,
                    ref_ligand_mol2=self.ref_ligand_mol2,
                    ddg=self.ddg,
                    rmsd=self.rmsd)

            if ddg_score:
                scores = ddg_score
            if len(smilies) == len(scores):
                pass
            else:
                scores = [10.0] * len(smilies)
        except:
            scores = [10.0] * len(smilies)
        
        self._internal_step += 1

        return ComponentResults([scores])

