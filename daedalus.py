"""
ARIADNE - DAEDALUS module

Functions for working with DAEDALUS outputs.
"""

import os
from typing import Tuple

from Bio.PDB.Model import Model
import Bio.SeqUtils.MeltingTemp as mt
import numpy as np
import pandas as pd

import cando
import pdb_tools
import seq_utils
from terms import *
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
        feasible_versions = [1, 2]
    else:
        feasible_versions = [version]
    # If no version is specified, infer the version from the txt file.
    succeeded_versions = list()
    failed_versions = list()
    for version in feasible_versions:
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
            stap_1xovers_scaf_nums_vis = base_annotations_vis.get((STAP_XO_1, rev), set())
            stap_1xovers_stap_nums_vis = cando.switch_strand(stap_1xovers_scaf_nums_vis, g_ax)
            stap_2xovers_scaf_nums_vis = base_annotations_vis.get((STAP_XO_2, rev), set())
            stap_2xovers_stap_nums_vis = cando.switch_strand(stap_2xovers_scaf_nums_vis, g_ax)
            all_xovers_scaf_nums_cando = (base_annotations_cando.get((SCAF, XO_1, rev, True), set()) |
                                          base_annotations_cando.get((SCAF, XO_1, fwd, False), set()) |
                                          base_annotations_cando.get((SCAF, XO_2, rev, True), set()) |
                                          base_annotations_cando.get((SCAF, XO_2, fwd, False), set()))
            all_xovers_stap_nums_cando = (base_annotations_cando.get((STAP, XO_1, fwd, True), set()) |
                                          base_annotations_cando.get((STAP, XO_1, rev, False), set()) |
                                          base_annotations_cando.get((STAP, XO_2, fwd, True), set()) |
                                          base_annotations_cando.get((STAP, XO_2, rev, False), set()))
            assert not (stap_1xovers_scaf_nums_vis | stap_2xovers_scaf_nums_vis) - all_xovers_scaf_nums_cando
            assert not (stap_1xovers_stap_nums_vis | stap_2xovers_stap_nums_vis) - all_xovers_stap_nums_cando
            base_annotations[SCAF, SCAF_XO, fwd] = all_xovers_scaf_nums_cando - stap_1xovers_scaf_nums_vis - stap_2xovers_scaf_nums_vis
            base_annotations[SCAF, STAP_XO_1, fwd] = stap_1xovers_scaf_nums_vis
            base_annotations[SCAF, STAP_XO_2, fwd] = stap_2xovers_scaf_nums_vis
            base_annotations[STAP, SCAF_XO, rev] = all_xovers_stap_nums_cando - stap_1xovers_stap_nums_vis - stap_2xovers_stap_nums_vis
            base_annotations[STAP, STAP_XO_1, rev] = stap_1xovers_stap_nums_vis
            base_annotations[STAP, STAP_XO_2, rev] = stap_2xovers_stap_nums_vis
        elif version == 2:
            scaf_xovers_scaf_nums_vis = base_annotations_vis.get((SCAF_XO, fwd), set())
            scaf_xovers_stap_nums_vis = cando.switch_strand(scaf_xovers_scaf_nums_vis, g_ax)
            all_1xovers_scaf_nums_cando = (base_annotations_cando.get((SCAF, XO_1, rev, True), set()) |
                                          base_annotations_cando.get((SCAF, XO_1, fwd, False), set()))
            all_1xovers_stap_nums_cando = (base_annotations_cando.get((STAP, XO_1, fwd, True), set()) |
                                          base_annotations_cando.get((STAP, XO_1, rev, False), set()))
            all_2xovers_scaf_nums_cando = (base_annotations_cando.get((SCAF, XO_2, rev, True), set()) |
                                          base_annotations_cando.get((SCAF, XO_2, fwd, False), set()))
            all_2xovers_stap_nums_cando = (base_annotations_cando.get((STAP, XO_2, fwd, True), set()) |
                                          base_annotations_cando.get((STAP, XO_2, rev, False), set()))                              
            assert not scaf_xovers_scaf_nums_vis - all_1xovers_scaf_nums_cando - all_2xovers_scaf_nums_cando
            assert not scaf_xovers_stap_nums_vis - all_1xovers_stap_nums_cando - all_2xovers_stap_nums_cando
            base_annotations[SCAF, SCAF_XO, fwd] = scaf_xovers_scaf_nums_vis
            base_annotations[SCAF, STAP_XO_1, fwd] = all_1xovers_scaf_nums_cando 
            base_annotations[SCAF, STAP_XO_2, fwd] = all_2xovers_scaf_nums_cando - scaf_xovers_scaf_nums_vis
            base_annotations[STAP, SCAF_XO, rev] = scaf_xovers_stap_nums_vis
            base_annotations[STAP, STAP_XO_1, rev] = all_1xovers_stap_nums_cando
            base_annotations[STAP, STAP_XO_2, rev] = all_2xovers_stap_nums_cando - scaf_xovers_stap_nums_vis
        else:
            raise ValueError(version)
        
        # bases adjacent to vertices (at edge ends)
        vertex_scaf_nums_vis = base_annotations_vis.get((VERTEX, fwd), set())
        vertex_stap_nums_vis = cando.switch_strand(vertex_scaf_nums_vis, g_ax)
        vertex_scaf_nums_cando = base_annotations_cando.get((SCAF, VERTEX, fwd, False), set())
        vertex_stap_nums_cando = base_annotations_cando.get((STAP, VERTEX, rev, False), set())
        assert vertex_scaf_nums_vis == vertex_scaf_nums_cando
        assert vertex_stap_nums_vis == vertex_stap_nums_cando
        base_annotations[SCAF, VERTEX, fwd] = vertex_scaf_nums_vis
        base_annotations[STAP, VERTEX, rev] = vertex_stap_nums_vis
    
    # strand termini
    scaf_term5_scaf_nums_vis = base_annotations_vis.get((SCAF_TM5, 0), set())
    scaf_term3_scaf_nums_vis = base_annotations_vis.get((SCAF_TM3, 0), set())
    scaf_term5_stap_nums_vis = cando.switch_strand(scaf_term5_scaf_nums_vis, g_ax)
    scaf_term3_stap_nums_vis = cando.switch_strand(scaf_term3_scaf_nums_vis, g_ax)
    stap_term5_scaf_nums_vis = base_annotations_vis.get((STAP_TM5, 0), set())
    stap_term3_scaf_nums_vis = base_annotations_vis.get((STAP_TM3, 0), set())
    stap_term5_stap_nums_vis = cando.switch_strand(stap_term5_scaf_nums_vis, g_ax)
    stap_term3_stap_nums_vis = cando.switch_strand(stap_term3_scaf_nums_vis, g_ax)
    scaf_term5_scaf_nums_cando = base_annotations_cando.get((SCAF, TM5, 0, False), set())
    scaf_term3_scaf_nums_cando = base_annotations_cando.get((SCAF, TM3, 0, False), set())
    scaf_term5_stap_nums_cando = base_annotations_cando.get((STAP, TM5, 0, True), set())
    scaf_term3_stap_nums_cando = base_annotations_cando.get((STAP, TM3, 0, True), set())
    stap_term5_scaf_nums_cando = base_annotations_cando.get((SCAF, TM5, 0, True), set())
    stap_term3_scaf_nums_cando = base_annotations_cando.get((SCAF, TM3, 0, True), set())
    stap_term5_stap_nums_cando = base_annotations_cando.get((STAP, TM5, 0, False), set())
    stap_term3_stap_nums_cando = base_annotations_cando.get((STAP, TM3, 0, False), set())
    assert scaf_term5_scaf_nums_vis == scaf_term5_scaf_nums_cando
    assert scaf_term3_scaf_nums_vis == scaf_term3_scaf_nums_cando
    assert scaf_term5_stap_nums_vis == scaf_term5_stap_nums_cando
    assert scaf_term3_stap_nums_vis == scaf_term3_stap_nums_cando
    base_annotations[SCAF, SCAF_TM5, 0] = scaf_term5_scaf_nums_vis
    base_annotations[SCAF, SCAF_TM3, 0] = scaf_term3_scaf_nums_vis
    base_annotations[SCAF, STAP_TM5, 0] = stap_term5_scaf_nums_cando
    base_annotations[SCAF, STAP_TM3, 0] = stap_term3_scaf_nums_cando
    base_annotations[STAP, SCAF_TM5, 0] = scaf_term5_stap_nums_vis
    base_annotations[STAP, SCAF_TM3, 0] = scaf_term3_stap_nums_vis
    base_annotations[STAP, STAP_TM5, 0] = stap_term5_stap_nums_cando
    base_annotations[STAP, STAP_TM3, 0] = stap_term3_stap_nums_cando
    
    # strand termini adjacent to single crossovers -- only defined from cando parser
    stap_term5XO_scaf_nums_cando = base_annotations_cando.get((SCAF, TM5_XO, 0, True), set())
    stap_term3XO_scaf_nums_cando = base_annotations_cando.get((SCAF, TM3_XO, 0, True), set())
    stap_term5XO_stap_nums_cando = base_annotations_cando.get((STAP, TM5_XO, 0, False), set())      
    stap_term3XO_stap_nums_cando = base_annotations_cando.get((STAP, TM3_XO, 0, False), set())      
    base_annotations[SCAF, STAP_TM5_XO, 0] = stap_term5XO_scaf_nums_cando
    base_annotations[SCAF, STAP_TM3_XO, 0] = stap_term3XO_scaf_nums_cando
    base_annotations[STAP, STAP_TM5_XO, 0] = stap_term5XO_stap_nums_cando
    base_annotations[STAP, STAP_TM3_XO, 0] = stap_term3XO_stap_nums_cando
    
    # mid-strand bases
    mid_scaf_nums_vis = base_annotations_vis.get(MIDDLE, set())
    mid_stap_nums_vis = cando.switch_strand(mid_scaf_nums_vis, g_ax)
    mid_scaf_nums_cando = base_annotations_cando.get((SCAF, MIDDLE, 0, False), set())
    mid_stap_nums_cando = base_annotations_cando.get((STAP, MIDDLE, 0, False), set())
    
    # the two parsers will not agree because the vis parser can't find staple crossovers and classifies them as middle
    base_annotations[SCAF, MIDDLE, 0] = mid_scaf_nums_cando
    base_annotations[STAP, MIDDLE, 0] = mid_stap_nums_cando
    
    # bases in vertices
    invertex_stap_nums_cando = base_annotations_cando.get((STAP, IN_VERTEX, 0, False), set())
    base_annotations[STAP, IN_VERTEX, 0] = invertex_stap_nums_cando
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


def walk_segment(base_info_df, g_up, g_dn, base_num_start):
    segment_5 = list(reversed(walk_segment_one_way(base_info_df, g_up, g_dn, "5", base_num_start)))
    segment_3 = walk_segment_one_way(base_info_df, g_up, g_dn, "3", base_num_start)
    assert base_num_start == segment_5[-1] == segment_3[0]
    segment = segment_5 + segment_3[1:]
    assert segment == sorted(segment) or segment == list(reversed(sorted(segment)))
    return segment


def get_next_base(base, g_up, g_dn, walk_direction):
    if walk_direction == "5":
        return g_up.get(base)
    elif walk_direction == "3":
        return g_dn.get(base)
    else:
        raise ValueError(walk_direction)


def walk_segment_one_way(base_info_df, g_up, g_dn, walk_direction, base_num_start):
    opp_direction = seq_utils.switch_direction(walk_direction)
    base_nums = [base_num_start]
    strand, feature, direction = base_info_df.loc[base_num_start, "Feature"].split("_")
    is_vertex = feature == IN_VERTEX
    while True:
        base = get_next_base(base_nums[-1], g_up, g_dn, walk_direction)
        if is_vertex:
            strand, feature, direction = base_info_df.loc[base, "Feature"].split("_")
            if feature != IN_VERTEX:
                return base_nums
        else:
            if feature in [VERTEX, SCAF_XO, STAP_XO_1, STAP_XO_2] and direction == opp_direction:
                return base_nums
            if feature in [STAP_TM5, STAP_TM5_XO, SCAF_TM3] and ((strand == STAP and walk_direction == "5") or (strand == SCAF and walk_direction == "3")):
                return base_nums
            if feature in [STAP_TM3, STAP_TM3_XO, SCAF_TM5] and ((strand == STAP and walk_direction == "3") or (strand == SCAF and walk_direction == "5")):
                return base_nums
            strand, feature, direction = base_info_df.loc[base, "Feature"].split("_")
        base_nums.append(base)


BOND_ID_VARS = ["Design", "CanDo number", "PDB chain", "PDB number", "Base", "Feature"]
BOND_TYPE_VARS = [f"{a1}-{a2} {dist}" for (n1, a1), (n2, a2) in pdb_tools.BACKBONE_BONDS for dist in ["length", "axial", "planar"]]
TM_SALTCORR = 6


def get_segments(base_info_df, g_up, g_dn):
    base_num = 0
    base_num_max = base_info_df.index.max()
    segments = list()
    while base_num < base_num_max:
        segment = walk_segment(base_info_df, g_up, g_dn, base_num + 1)
        seg5 = segment[0]
        seg3 = segment[-1]
        segments.append({
            "Seg5": seg5,
            "Seg3": seg3,
            "Feature5": base_info_df.loc[seg5, "Feature"],
            "Feature3": base_info_df.loc[seg3, "Feature"],
        })
        base_num = seg3
    segments = pd.DataFrame(segments)
    return segments


def get_melt(seq):
    return mt.Tm_NN(seq, nn_table=mt.R_DNA_NN1, saltcorr=6)


def assemble_base_info(design, model, base_annotations, cando_num_to_pdb_chain_num, pair_directions, g_up, g_dn, g_ax, base_seq, compute_bond_lengths=True):
    base_num_to_annotation = {base_num: annotation for annotation, base_nums in base_annotations.items() for base_num in
                              base_nums}
    pdb_chain_num_to_cando_num = {pdb_chain_num: base_num for base_num, pdb_chain_num in
                                  cando_num_to_pdb_chain_num.items()}
    if compute_bond_lengths:
        bond_types_and_lengths, bond_types_and_axial_distances, bond_types_and_planar_distances = pdb_tools.get_bond_types_and_lengths(model, pdb_chain_num_to_cando_num, pair_directions, base_num_to_annotation, g_ax)
        bond_info_fields = BOND_ID_VARS + BOND_TYPE_VARS
    else:
        bond_info_fields = BOND_ID_VARS
    bases_info = list()
    for chain in model:
        pdb_chain = chain.get_id()
        for residue in chain:
            pdb_num = residue.get_id()[1]
            base_num = pdb_chain_num_to_cando_num[pdb_chain, pdb_num]
            strand, feature, direction = base_num_to_annotation[base_num]
            annotation_label = f"{strand}_{feature}_{direction}"
            # Add the basic base annotations.
            base_info = {"Design": design,
                         "CanDo number": base_num,
                         "PDB chain": pdb_chain,
                         "PDB number": pdb_num,
                         "Base": base_seq[base_num],
                         "Feature": annotation_label}
            if compute_bond_lengths:
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
    base_info_df = pd.DataFrame.from_records(bases_info, columns=bond_info_fields)
    base_info_df.index = base_info_df["CanDo number"]
    # Compute segment-based properties
    segment_5ps = list()
    segment_3ps = list()
    segment_seqs = list()
    features_end_5 = list()
    features_end_3 = list()
    for base_num in base_info_df.index:
        segment = walk_segment(base_info_df, g_up, g_dn, base_num)
        segment_5ps.append(segment[0])
        segment_3ps.append(segment[-1])
        segment_seq = "".join([base_info_df.loc[base_num, "Base"] for base_num in segment])
        segment_seqs.append(segment_seq)
        strand, feature_end_5, direction = base_info_df.loc[segment[0], "Feature"].split("_")
        features_end_5.append(feature_end_5)
        strand, feature_end_3, direction = base_info_df.loc[segment[-1], "Feature"].split("_")
        features_end_3.append(feature_end_3)
        seg_features = [base_info_df.loc[x, "Feature"].split("_")[1] for x in segment]
        non_middle = [f for f in seg_features[1: -1] if f != "Middle"]
    base_info_df["Segment5'"] = segment_5ps
    base_info_df["Segment3'"] = segment_3ps
    #base_info_df["SegmentSeq"] = segment_seqs
    base_info_df["SegmentLength"] = list(map(len, segment_seqs))
    assert np.all(base_info_df["SegmentLength"] == np.abs(base_info_df["Segment3'"] - base_info_df["Segment5'"]) + 1)
    #base_info_df["SegmentGC"] = list(map(seq_utils.get_gc_content, segment_seqs))
    #base_info_df["SegmentTm"] = list(map(get_melt, segment_seqs))
    base_info_df["Segment5'Feature"] = features_end_5
    base_info_df["Segment3'Feature"] = features_end_3
    return base_info_df


def get_unique_segments(base_info_df):
    unique_segments = sorted({(base_info_df.loc[base_num, "Segment5'"], base_info_df.loc[base_num, "Segment3'"]) for base_num in base_info_df.index})
    return unique_segments


def get_unique_segment_seqs(base_info_df):
    bases = base_info_df["Base"]
    unique_segments = get_unique_segments(base_info_df)
    unique_segment_seqs = {(seg_5p, seg_3p): "".join(bases.loc[seg_5p: seg_3p]) for seg_5p, seg_3p in unique_segments}
    return unique_segment_seqs


def get_segment_seqs(base_info_df, update=False):
    unique_segment_seqs = get_unique_segment_seqs(base_info_df)
    all_segment_seqs = pd.Series([unique_segment_seqs[seg_5p, seg_3p] for seg_5p, seg_3p in zip(base_info_df["Segment5'"], base_info_df["Segment3'"])], index=base_info_df.index)
    if update:
        base_info_df["SegmentSeq"] = all_segment_seqs
    return all_segment_seqs


def get_segment_end_dists(base_info_df):
    base_nums = base_info_df.index
    dist_5p = base_nums - base_info_df["Segment5'"]
    dist_3p = base_info_df["Segment3'"] - base_nums
    assert all(dist_5p >= 0)
    assert all(dist_3p >= 0)
    return dist_5p, dist_3p


def get_weighted_gc(base_info_df, weight_func, update=False):
    weighted_gcs = pd.Series(index=base_info_df.index)
    weighted_gcs.name = "WeightedGC"
    segments_seqs = get_unique_segment_seqs(base_info_df)
    segments_base_nums = {(seg_5p, seg_3p): np.array(list(range(seg_5p, seg_3p + 1))) for seg_5p, seg_3p in segments_seqs}
    segments_is_gc = {(seg_5p, seg_3p): np.array([base in "GC" for base in seq]) for (seg_5p, seg_3p), seq in segments_seqs.items()}
    for base_num in base_info_df.index:
        seg_5p = base_info_df.loc[base_num, "Segment5'"]
        seg_3p = base_info_df.loc[base_num, "Segment3'"]
        is_gc = segments_is_gc[seg_5p, seg_3p]
        seg_nums = segments_base_nums[seg_5p, seg_3p]
        weights = list(map(weight_func, seg_nums - base_num))
        weight_sum = sum(weights)
        if weight_sum <= 0:
            raise ValueError("weight_sum must be positive.")
        weighted_gc = np.dot(weights, is_gc) / weight_sum
        weighted_gcs[base_num] = weighted_gc
    return weighted_gcs


def get_feature_type(full_feature):
    strand, feature_type, direction = full_feature.split("_")
    return feature_type


adirectional_features = [MIDDLE]
unidirectional_features = {"5": [STAP_TM5, STAP_TM5_XO, SCAF_TM5],
                           "3": [STAP_TM3, STAP_TM3_XO, SCAF_TM3]}
bidirectional_features = [SCAF_XO, STAP_XO_1, STAP_XO_2, VERTEX]

def format_feature_name(feature, direction):
    direction = str(direction)
    if direction in ["5", "3"]:
        direction = {"5": "before", "3": "after"}[direction]
        return f"{direction}_{feature}"
    elif direction in ["0"]:
        if feature in bidirectional_features:
            raise ValueError(f"{feature} with direction {direction}")
        for direction in ["5", "3"]:
            if feature in unidirectional_features[direction]:
                return format_feature_name(feature, direction)
        if feature in adirectional_features:
            return feature
        else:
            raise ValueError(feature)
    else:
        raise ValueError(direction)


features_order = [SCAF_XO, STAP_XO_2, STAP_XO_1, VERTEX, STAP_TM5, STAP_TM3, STAP_TM5_XO, STAP_TM3_XO]

def get_formatted_features(base_info_df):
    base_info_features = set(base_info_df["Feature"])
    formatted_features = list()
    for feature in features_order:
        if feature in bidirectional_features:
            directions = ["5", "3"]
        elif feature in unidirectional_features["5"]:
            directions = ["5"]
        elif feature in unidirectional_features["3"]:
            directions = ["3"]
        elif feature in adirectional_features:
            continue
        else:
            raise ValueError("{feature} not found")
        for direction in directions:
            feature_label = f"{SCAF}_{feature}"
            if any((x.startswith(feature_label) for x in base_info_features)):
                formatted_features.append(format_feature_name(feature, direction))
    return formatted_features


def get_coeffs_matrix(base_info_df, funcs):
    dist_5p, dist_3p = get_segment_end_dists(base_info_df)
    features = get_formatted_features(base_info_df)
    coeffs = pd.DataFrame(index=base_info_df.index, columns=features)
    for base_num in base_info_df.index:
        # Here the direction refers to the direction from the feature to the base.
        # The base of interest is in the 3' direction from the 5' feature,
        # so the direction is 3.
        feature_5p = format_feature_name(base_info_df.loc[base_num, "Segment5'Feature"], "3")
        if feature_5p in features:
            coeff_5p = funcs[feature_5p](dist_5p.loc[base_num])
            coeffs.loc[base_num, feature_5p] = coeff_5p
        # The base of interest is in the 5' direction from the 3' feature,
        # so the direction is 5.
        feature_3p = format_feature_name(base_info_df.loc[base_num, "Segment3'Feature"], "5")
        if feature_3p in features:
            coeff_3p = funcs[feature_3p](dist_3p.loc[base_num])
            coeffs.loc[base_num, feature_3p] = coeff_3p
    coeffs.fillna(0, inplace=True)
    return coeffs

