"""
ARIADNE - PDB module

Functions for working with PDBs and their Biopython representations.
"""

from collections import defaultdict
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd

from Bio.PDB import PDBParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure


BACKBONE_BONDS = [
        ((-1, "O3'"), (0, "P")),
        ((0, "P"), (0, "O5'")),
        ((0, "O5'"), (0, "C5'")),
        ((0, "C5'"), (0, "C4'")),
        ((0, "C4'"), (0, "C3'")),
        ((0, "C3'"), (0, "O3'")),
]


def get_structure(pdb_file, name="structure"):
    structure = PDBParser().get_structure(name, pdb_file)
    return structure


def get_model(structure: Union[Model, Structure]) -> Model:
    """
    Get the first model from a Biopython PDB structure.

    :param structure: structure or model object
    :return model: the first model in the structure
    """
    if isinstance(structure, Model):
        model = structure
    else:
        models = list(structure.get_models())
        model = models[0]
    return model


def get_chain_nums_seq(chain: Chain, na: Optional[str] = None) -> Tuple[List[int], str]:
    """
    Get the sequence and base numbers of a chain.

    :param chain: Biopython chain from which to extract the sequence
           na: whether to convert the sequence to DNA or RNA (optional)
    :return nums: the number of each residue in the chain
            seq: the sequence of the chain using one-letter codes
    """
    nums = [residue.get_id()[1] for residue in chain]
    seq = "".join([residue.get_resname()[0] for residue in chain])
    if na is not None:
        if na.upper() == "DNA":
            seq = seq.replace("U", "T")
        elif na.upper() == "RNA":
            seq = seq.replace("T", "U")
        else:
            raise ValueError(na)
    return nums, seq


def get_chains_nums_seqs(chains: List[Chain], na: Optional[str] = None) -> Dict[str, Tuple[List[int], str]]:
    """
    Get the sequence and base numbers of every chain in a structure or model.
    :param chains:
    :return:
    """
    chains_nums_seqs = {chain.get_id(): get_chain_nums_seq(chain, na) for chain in chains}
    return chains_nums_seqs


def get_bond_types_and_lengths(model):
    bond_types_and_lengths = dict()
    for chain in model:
        chain_id = chain.get_id()
        prev_base = None
        for base in chain:
            base_num = base.get_id()[1]
            bond_types_and_lengths[(chain_id, base_num)] = dict()
            for bond_type in BACKBONE_BONDS:
                (num1_del, atom1_name), (num2_del, atom2_name) = bond_type
                if num1_del == 0:
                    atom1 = base[atom1_name]
                elif num1_del == -1:
                    atom1 = prev_base[atom1_name] if prev_base is not None else None
                else:
                    raise ValueError(num1_del)
                if num2_del == 0:
                    atom2 = base[atom2_name]
                elif num2_del == -1:
                    atom2 = prev_base[atom2_name] if prev_base is not None else None
                else:
                    raise ValueError(num1_del)
                if atom1 and atom2:
                    length = atom2 - atom1
                    bond_types_and_lengths[(chain_id, base_num)][bond_type] = length
            prev_base = base
    return bond_types_and_lengths
