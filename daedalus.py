"""
ARIADNE - DAEDALUS module

Functions for working with DAEDALUS outputs.
"""

import os
from typing import Tuple

from Bio.PDB.Model import Model
import Bio.SeqUtils.MeltingTemp as mt
import pandas as pd

import cando
import pdb_tools
import seq_utils
import terms
import vis


def get_files(design_directory: str) -> Tuple[str, str, str]:
    """
    Locate and return three files in the directory containing the outputs of a DAEDALUS design.

    :param design_directory: path to the directory of the design
    :return pdb_file: the path to the multi-chain PDB file
    :return cando_file: the path to the CanDo file
    :return vis_file: the path to the ASCII visualization file
    """
    pdb_file, vis_file, cando_file = None, None, None
    for fname in os.listdir(design_directory):
        fpath = os.path.join(design_directory, fname)
        if os.path.splitext(fname)[1] == ".pdb" and "chseg" not in fname and "multimodel" not in fname:
            assert pdb_file is None
            pdb_file = fpath
        elif os.path.splitext(fname)[1] == ".cndo":
            assert cando_file is None
            cando_file = fpath
        elif os.path.splitext(fname)[1] == ".txt" and fname.startswith("seq_"):
            assert vis_file is None
            vis_file = fpath
    assert pdb_file is not None and vis_file is not None and cando_file is not None
    return pdb_file, cando_file, vis_file


def annotate_bases(g_up, g_dn, g_ax, vis_file, version=None):
    """
    Find base annotations from both CanDo and ASCII visualization files.
    :param g_up:
    :param g_dn:
    :param g_ax:
    :param vis_file:
    :return:
    """
    if version is None:
        # If no version is specified, infer the version from the txt file.
        succeeded_versions = list()
        failed_versions = list()
        for version in [1, 2]:
            try:
                base_annotations_vis, edges = vis.annotate_bases(vis_file, version)
            except AssertionError:
                failed_versions.append(version)
            else:
                succeeded_versions.append(version)
        assert len(succeeded_versions) == 1
        version = succeeded_versions[0]
    base_nums = sorted(g_up)
    assert base_nums == sorted(g_dn) == sorted(g_ax)
    # Get the base annotations from CanDo and the ASCII vis file.
    base_annotations_cando = cando.annotate_bases(g_up, g_dn, g_ax)
    # Refine the annotation for each base.
    base_annotations = dict()
    for fwd, rev in zip([5, 3], [3, 5]):
        # crossovers
        if version == 1:
            stap_xovers_scaf_nums_vis = base_annotations_vis.get((terms.STAP_XO, rev), set())
            stap_xovers_stap_nums_vis = cando.switch_strand(stap_xovers_scaf_nums_vis, g_ax)
            all_xovers_scaf_nums_cando = (base_annotations_cando.get((terms.SCAF, terms.XO, rev, True), set()) |
                                          base_annotations_cando.get((terms.SCAF, terms.XO, fwd, False), set()))
            all_xovers_stap_nums_cando = (base_annotations_cando.get((terms.STAP, terms.XO, fwd, True), set()) |
                                          base_annotations_cando.get((terms.STAP, terms.XO, rev, False), set()))
            assert not stap_xovers_scaf_nums_vis - all_xovers_scaf_nums_cando
            assert not stap_xovers_stap_nums_vis - all_xovers_stap_nums_cando
            base_annotations[terms.SCAF, terms.SCAF_XO, fwd] = all_xovers_scaf_nums_cando - stap_xovers_scaf_nums_vis
            base_annotations[terms.SCAF, terms.STAP_XO, fwd] = stap_xovers_scaf_nums_vis
            base_annotations[terms.STAP, terms.SCAF_XO, rev] = all_xovers_stap_nums_cando - stap_xovers_stap_nums_vis
            base_annotations[terms.STAP, terms.STAP_XO, rev] = stap_xovers_stap_nums_vis
        elif version == 2:
            scaf_xovers_scaf_nums_vis = base_annotations_vis.get((terms.SCAF_XO, fwd), set())
            scaf_xovers_stap_nums_vis = cando.switch_strand(scaf_xovers_scaf_nums_vis, g_ax)
            all_xovers_scaf_nums_cando = (base_annotations_cando.get((terms.SCAF, terms.XO, rev, True), set()) |
                                          base_annotations_cando.get((terms.SCAF, terms.XO, fwd, False), set()))
            all_xovers_stap_nums_cando = (base_annotations_cando.get((terms.STAP, terms.XO, fwd, True), set()) |
                                          base_annotations_cando.get((terms.STAP, terms.XO, rev, False), set()))
            assert not scaf_xovers_scaf_nums_vis - all_xovers_scaf_nums_cando
            assert not scaf_xovers_stap_nums_vis - all_xovers_stap_nums_cando
            base_annotations[terms.SCAF, terms.SCAF_XO, fwd] = scaf_xovers_scaf_nums_vis
            base_annotations[terms.SCAF, terms.STAP_XO, fwd] = all_xovers_scaf_nums_cando - scaf_xovers_scaf_nums_vis
            base_annotations[terms.STAP, terms.SCAF_XO, rev] = scaf_xovers_stap_nums_vis
            base_annotations[terms.STAP, terms.STAP_XO, rev] = all_xovers_stap_nums_cando - scaf_xovers_stap_nums_vis
        else:
            raise ValueError(version)
        # strand termini
        scaf_term_scaf_nums_vis = base_annotations_vis.get((terms.SCAF_TM, fwd), set())
        scaf_term_stap_nums_vis = cando.switch_strand(scaf_term_scaf_nums_vis, g_ax)
        stap_term_scaf_nums_vis = base_annotations_vis.get((terms.STAP_TM, fwd), set())
        stap_term_stap_nums_vis = cando.switch_strand(stap_term_scaf_nums_vis, g_ax)
        scaf_term_scaf_nums_cando = base_annotations_cando.get((terms.SCAF, terms.TM, fwd, False), set())
        scaf_term_stap_nums_cando = base_annotations_cando.get((terms.STAP, terms.TM, fwd, True), set())
        stap_term_scaf_nums_cando = base_annotations_cando.get((terms.SCAF, terms.TM, fwd, True), set())
        stap_term_stap_nums_cando = base_annotations_cando.get((terms.STAP, terms.TM, fwd, False), set())
        assert scaf_term_scaf_nums_vis == scaf_term_scaf_nums_cando
        assert scaf_term_stap_nums_vis == scaf_term_stap_nums_cando
        assert stap_term_scaf_nums_vis == stap_term_scaf_nums_cando
        assert stap_term_stap_nums_vis == stap_term_stap_nums_cando
        base_annotations[terms.SCAF, terms.SCAF_TM, fwd] = scaf_term_scaf_nums_vis
        base_annotations[terms.SCAF, terms.STAP_TM, fwd] = stap_term_scaf_nums_vis
        base_annotations[terms.STAP, terms.SCAF_TM, fwd] = scaf_term_stap_nums_vis
        base_annotations[terms.STAP, terms.STAP_TM, fwd] = stap_term_stap_nums_vis
        # edge ends
        edge_end_scaf_nums_vis = base_annotations_vis.get((terms.EDGE_TM, fwd), set())
        edge_end_stap_nums_vis = cando.switch_strand(edge_end_scaf_nums_vis, g_ax)
        edge_end_scaf_nums_cando = base_annotations_cando.get((terms.SCAF, terms.EDGE_TM, fwd, False), set())
        edge_end_stap_nums_cando = base_annotations_cando.get((terms.STAP, terms.EDGE_TM, rev, False), set())
        assert edge_end_scaf_nums_vis == edge_end_scaf_nums_cando
        assert edge_end_stap_nums_vis == edge_end_stap_nums_cando
        base_annotations[terms.SCAF, terms.EDGE_TM, fwd] = edge_end_scaf_nums_vis
        base_annotations[terms.STAP, terms.EDGE_TM, rev] = edge_end_stap_nums_vis
    # mid-strand bases
    mid_scaf_nums_vis = base_annotations_vis.get(terms.MIDDLE, set())
    mid_stap_nums_vis = cando.switch_strand(mid_scaf_nums_vis, g_ax)
    mid_scaf_nums_cando = base_annotations_cando.get((terms.SCAF, terms.MIDDLE, 0, False), set())
    mid_stap_nums_cando = base_annotations_cando.get((terms.STAP, terms.MIDDLE, 0, False), set())
    # the two parsers will not agree because the vis parser can't find staple crossovers and classifies them as middle
    base_annotations[terms.SCAF, terms.MIDDLE, 0] = mid_scaf_nums_cando
    base_annotations[terms.STAP, terms.MIDDLE, 0] = mid_stap_nums_cando
    # vertex bases
    vert_stap_nums_cando = base_annotations_cando.get((terms.STAP, terms.VERTEX, 0, False), set())
    base_annotations[terms.STAP, terms.VERTEX, 0] = vert_stap_nums_cando
    base_annotation_nums = sorted([base_num for annot_nums in base_annotations.values() for base_num in annot_nums])
    assert base_annotation_nums == base_nums
    return base_annotations, edges, version


def map_cando_num_to_pdb_chain_num(model: Model, base_annotations, g_dn, g_ax):
    """
    :return:
    """
    # Get sequences and base numbers from PDB.
    scaffold_chain = model[seq_utils.scaffold_chain_id]
    scaffold_base_nums_pdb, scaffold_seq_pdb = pdb_tools.get_chain_nums_seq(scaffold_chain, na="RNA")
    staples_chains_nums_seqs_pdb = pdb_tools.get_chains_nums_seqs(
        [chain for chain in model if chain.get_id() != seq_utils.scaffold_chain_id], na="DNA")
    staples_seqs_pdb = {seq: (chain, nums) for chain, (nums, seq) in staples_chains_nums_seqs_pdb.items()}
    # Get staple sequences and base numbers from CanDo.
    staples_bases_nums_cando = cando.get_staples_bases_nums(base_annotations, g_dn)
    staples_seqs_cando = cando.get_staples_seqs(staples_bases_nums_cando, g_ax, scaffold_seq_pdb)
    # Map the CanDo file index to the pdb file chain and number.
    assert sorted(staples_seqs_pdb) == sorted(staples_seqs_cando)
    cando_num_to_pdb_chain_num = {base_num: (seq_utils.scaffold_chain_id, base_num) for base_num, base_seq in
                                  zip(scaffold_base_nums_pdb, scaffold_seq_pdb)}
    for staple_bases_nums_cando, staple_seq_cando in zip(staples_bases_nums_cando, staples_seqs_cando):
        chain_pdb, staple_bases_nums_pdb = staples_seqs_pdb[staple_seq_cando]
        for base_num_cando, base_num_pdb in zip(staple_bases_nums_cando, staple_bases_nums_pdb):
            assert base_num_cando not in cando_num_to_pdb_chain_num
            cando_num_to_pdb_chain_num[base_num_cando] = chain_pdb, base_num_pdb
    assert sorted(cando_num_to_pdb_chain_num) == list(range(min(cando_num_to_pdb_chain_num), max(cando_num_to_pdb_chain_num) + 1))
    return cando_num_to_pdb_chain_num


def walk_segment(bases_info_df, g_up, g_dn, base_num_start):
    segment_5 = list(reversed(walk_segment_one_way(bases_info_df, g_up, g_dn, "5", base_num_start)))
    segment_3 = walk_segment_one_way(bases_info_df, g_up, g_dn, "3", base_num_start)
    assert base_num_start == segment_5[-1] == segment_3[0]
    segment = segment_5 + segment_3[1:]
    return segment


def walk_segment_one_way(bases_info_df, g_up, g_dn, walk_direction, base_num_start):
    base_nums = list()
    current_base = base_num_start
    strand, feature, direction = bases_info_df.loc[current_base, "location"].split("_")
    vertex = feature == terms.VERTEX
    is_segment_end = False
    while not is_segment_end:
        base_nums.append(current_base)
        if not vertex:
            is_scaf_terminus = (strand == terms.SCAF and (
                    (terms.SCAF_TM == feature and direction == walk_direction) or
                    (terms.STAP_TM == feature and direction == seq_utils.switch_direction(walk_direction))))
            is_stap_terminus = (strand == terms.STAP and (
                    (terms.STAP_TM == feature and direction == walk_direction) or
                    (terms.SCAF_TM == feature and direction == seq_utils.switch_direction(walk_direction))))
            is_edge_terminus = terms.EDGE_TM == feature and direction == walk_direction
            is_crossover = terms.XO in feature and direction == seq_utils.switch_direction(walk_direction)
            is_segment_end = is_scaf_terminus or is_stap_terminus or is_edge_terminus or is_crossover
        if not is_segment_end:
            if walk_direction == "5":
                current_base = g_up.get(current_base)
            elif walk_direction == "3":
                current_base = g_dn.get(current_base)
            else:
                raise ValueError(walk_direction)
            strand, feature, direction = bases_info_df.loc[current_base, "location"].split("_")
        if vertex:
            is_segment_end = feature != terms.VERTEX
    return base_nums


BOND_ID_VARS = ["design", "CanDo number", "PDB chain", "PDB number", "base", "location"]
BOND_TYPE_VARS = [f"{a1}-{a2} {dist}" for (n1, a1), (n2, a2) in pdb_tools.BACKBONE_BONDS for dist in ["length", "axial", "planar"]]
BOND_INFO_FIELDS = BOND_ID_VARS + BOND_TYPE_VARS
TM_SALTCORR = 6


def get_melt(seq):
    return mt.Tm_NN(seq, nn_table=mt.R_DNA_NN1, saltcorr=6)


def assemble_base_info(design, model, base_annotations, cando_num_to_pdb_chain_num, pair_directions, g_up, g_dn, g_ax, base_seq):
    base_num_to_annotation = {base_num: annotation for annotation, base_nums in base_annotations.items() for base_num in
                              base_nums}
    pdb_chain_num_to_cando_num = {pdb_chain_num: base_num for base_num, pdb_chain_num in
                                  cando_num_to_pdb_chain_num.items()}
    bond_types_and_lengths, bond_types_and_axial_distances, bond_types_and_planar_distances = pdb_tools.get_bond_types_and_lengths(model, pdb_chain_num_to_cando_num, pair_directions, base_num_to_annotation, g_ax)
    bases_info = list()
    for (pdb_chain, pdb_num) in bond_types_and_lengths:
        base_num = pdb_chain_num_to_cando_num[pdb_chain, pdb_num]
        strand, feature, direction = base_num_to_annotation[base_num]
        annotation_label = f"{strand}_{feature}_{direction}"
        # Add the basic base annotations.
        base_info = {"design": design,
                     "CanDo number": base_num,
                     "PDB chain": pdb_chain,
                     "PDB number": pdb_num,
                     "base": base_seq[base_num],
                     "location": annotation_label}
        # Add the bond types and lengths
        base_bond_types_and_lengths = bond_types_and_lengths[pdb_chain, pdb_num]
        for bond_type, bond_length in base_bond_types_and_lengths.items():
            ((num1_del, atom1), (num2_del, atom2)) = bond_type
            bond_type_label = f"{atom1}-{atom2} length"
            base_info[bond_type_label] = bond_length
        base_bond_types_and_axial_distances = bond_types_and_axial_distances[pdb_chain, pdb_num]
        for bond_type, axial_dist in base_bond_types_and_axial_distances.items():
            ((num1_del, atom1), (num2_del, atom2)) = bond_type
            bond_type_label = f"{atom1}-{atom2} axial"
            base_info[bond_type_label] = axial_dist
        base_bond_types_and_planar_distances = bond_types_and_planar_distances[pdb_chain, pdb_num]
        for bond_type, planar_dist in base_bond_types_and_planar_distances.items():
            ((num1_del, atom1), (num2_del, atom2)) = bond_type
            bond_type_label = f"{atom1}-{atom2} planar"
            base_info[bond_type_label] = planar_dist
        bases_info.append(base_info)
    # Convert to dataframe
    bases_info_df = pd.DataFrame.from_records(bases_info, columns=BOND_INFO_FIELDS)
    bases_info_df.index = bases_info_df["CanDo number"]
    # Compute segment-based properties
    segment_seqs = list()
    features_end_5 = list()
    features_end_3 = list()
    for base_num in bases_info_df.index:
        segment = walk_segment(bases_info_df, g_up, g_dn, base_num)
        segment_seq = "".join([bases_info_df.loc[base_num, "base"] for base_num in segment])
        segment_seqs.append(segment_seq)
        strand, feature_end_5, direction = bases_info_df.loc[segment[0], "location"].split("_")
        features_end_5.append(feature_end_5)
        strand, feature_end_3, direction = bases_info_df.loc[segment[-1], "location"].split("_")
        features_end_3.append(feature_end_3)
    bases_info_df["SegmentSeq"] = segment_seqs
    bases_info_df["SegmentLength"] = list(map(len, segment_seqs))
    bases_info_df["SegmentGC"] = list(map(seq_utils.get_gc_content, segment_seqs))
    bases_info_df["SegmentMelt"] = list(map(get_melt, segment_seqs))
    bases_info_df["Segment5Feature"] = features_end_5
    bases_info_df["Segment3Feature"] = features_end_3
    return bases_info_df
