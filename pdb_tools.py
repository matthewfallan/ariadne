"""
ARIADNE - PDB module

Functions for working with PDBs and their Biopython representations.
"""

from collections import defaultdict
from typing import Dict, List, Optional, Tuple, Union

from Bio.PDB import PDBParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
import numpy as np
import pandas as pd

import terms


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


def get_bond_types_and_lengths(model, pdb_chain_num_to_cando_num, pair_directions, base_num_to_annotation, g_ax):
    bond_types_and_lengths = dict()
    bond_types_and_axial_distances = dict()
    bond_types_and_planar_distances = dict()
    for chain in model:
        chain_id = chain.get_id()
        prev_base = None
        for base in chain:
            base_num = base.get_id()[1]
            pdb_chain_num = chain_id, base_num
            cando_num = pdb_chain_num_to_cando_num[pdb_chain_num]
            paired_num = g_ax.get(cando_num)
            scaf_num = min(cando_num, paired_num) if paired_num is not None else None
            bond_types_and_lengths[pdb_chain_num] = dict()
            bond_types_and_axial_distances[pdb_chain_num] = dict()
            bond_types_and_planar_distances[pdb_chain_num] = dict()
            for bond_type in BACKBONE_BONDS:
                (num1_del, atom1_name), (num2_del, atom2_name) = bond_type
                if prev_base:
                    cando_prev = pdb_chain_num_to_cando_num[chain_id, prev_base.get_id()[1]]
                    paired_prev = g_ax.get(cando_prev)
                    scaf_prev = min(cando_prev, paired_prev) if paired_prev is not None else None
                else:
                    scaf_prev = None
                if num1_del == 0:
                    atom1 = base[atom1_name]
                    axis1 = pair_directions.loc[scaf_num] if scaf_num is not None else None
                elif num1_del == -1:
                    atom1 = prev_base[atom1_name] if prev_base is not None else None
                    axis1 = pair_directions.loc[scaf_prev] if scaf_prev is not None else None
                else:
                    raise ValueError(num1_del)
                if num2_del == 0:
                    atom2 = base[atom2_name]
                    axis2 = pair_directions.loc[scaf_num] if scaf_num is not None else None
                elif num2_del == -1:
                    atom2 = prev_base[atom2_name] if prev_base is not None else None
                    axis2 = pair_directions.loc[scaf_prev] if scaf_prev is not None else None
                else:
                    raise ValueError(num1_del)
                if atom1 and atom2:
                    # find the vector connecting the two atoms
                    atom_displacement = atom2.get_coord() - atom1.get_coord()
                    distance = np.linalg.norm(atom_displacement)
                    bond_types_and_lengths[pdb_chain_num][bond_type] = distance
                    if axis1 is None or axis2 is None:
                        axial_dist = np.nan
                    elif base_num_to_annotation[cando_num][1] == terms.EDGE_TM and (num1_del == -1 or num2_del == -1):
                        # axis switches at edge termini, so can't compute axial distance if using the previous base
                        axial_dist = np.nan
                    else:
                        # project the vector onto the axis of each helix
                        axial_disp_1 = project_vector(atom_displacement, axis1)
                        axial_dist_1 = np.linalg.norm(axial_disp_1)
                        axial_disp_2 = project_vector(atom_displacement, axis2)
                        axial_dist_2 = np.linalg.norm(axial_disp_2)
                        assert np.isclose(axial_dist_1, axial_dist_2)
                        axial_dist = (axial_dist_1 + axial_dist_2) / 2
                    # compute the planar distance using the Pythagorean theorem
                    # distance^2 = axial_dist^2 + planar_dist^2
                    planar_dist = np.sqrt(distance**2 - axial_dist**2)
                    bond_types_and_axial_distances[pdb_chain_num][bond_type] = axial_dist
                    bond_types_and_planar_distances[pdb_chain_num][bond_type] = planar_dist
            prev_base = base
    return bond_types_and_lengths, bond_types_and_axial_distances, bond_types_and_planar_distances


def project_vector(query, target):
    """
    Project vector query onto vector target.
    :param query:
    :param target:
    :param query_center:
    :param target_center:
    :return:
    """
    assert query.shape == target.shape
    projection = np.dot(query, target) / np.linalg.norm(target)
    return projection
