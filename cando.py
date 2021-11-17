"""
ARIADNE - CanDo Module

Functions for processing data in CanDo files.
"""

from collections import defaultdict
import itertools
from typing import Dict, List, Optional, Set, Tuple, Union

import pandas as pd

import seq_utils
from terms import *


def get_connectivity(cando_file: str) -> Tuple[
        List[int], Dict[int, Optional[int]],
        Dict[int, Optional[int]], Dict[int, Optional[int]]]:
    """
    Retrieve residue connectivity information from a CanDo file.

    :param cando_file: str, file path of CanDo file
    :return
    Returns
    -------
    base_nums : list[int]
        number of every base
    g_up : dict[int, (int, None)]
        the number of the upstream residue (if none, then None)
    g_dn : dict[int, (int, None)]
        the number of the downstream residue (if none, then None)
    g_ax : dict[int, (int, None)]
        the number of the paired residue (if none, then None)
    """
    base_seq = dict()
    g_up = dict()
    g_dn = dict()
    g_ax = dict()
    with open(cando_file) as f:
        # Read lines until reaching the header.
        header = "dnaTop,id,up,down,across,seq"
        line = f.readline()
        while line.strip() != header:
            line = f.readline()
        # Read the connectivity information.
        line = f.readline()
        while line.strip():
            # Read the information, which is comma-delimited.
            dna_top, num, up, down, across, seq = line.strip().split(",")
            dna_top, num, up, down, across = int(dna_top), int(num), int(up), int(down), int(across)
            # Add the information to the graphs.
            base_seq[num] = seq
            g_up[num] = up if up != -1 else None
            g_dn[num] = down if down != -1 else None
            g_ax[num] = across if across != -1 else None
            line = f.readline()
    # List all of the base numbers in ascending order.
    base_nums = sorted(g_up)
    return base_nums, base_seq, g_up, g_dn, g_ax


def get_pair_centers(cando_file: str):
    """
    Return the coordinates from a CanDo file.
    :param cando_file:
    :return:
    """
    xs = dict()
    ys = dict()
    zs = dict()
    with open(cando_file) as f:
        # Read lines until reaching the header.
        header = 'dNode,"e0(1)","e0(2)","e0(3)"'
        line = f.readline()
        while line.strip() != header:
            line = f.readline()
        # Read the connectivity information.
        line = f.readline()
        while line.strip():
            # Read the information, which is comma-delimited.
            num, x, y, z = line.strip().split(",")
            num, x, y, z = int(num), float(x), float(y), float(z)
            # Record the information.
            xs[num] = x
            ys[num] = y
            zs[num] = z
            line = f.readline()
    assert sorted(xs) == list(range(min(xs), max(xs) + 1))
    centers = pd.DataFrame.from_dict({"x": xs, "y": ys, "z": zs})
    return centers


def get_pair_directions(cando_file: str):
    """
    Return the directions of the pairs from a CanDo file.
    :param cando_file:
    :return:
    """
    xs = dict()
    ys = dict()
    zs = dict()
    with open(cando_file) as f:
        # Read lines until reaching the header.
        header = 'triad,"e1(1)","e1(2)","e1(3)","e2(1)","e2(2)","e2(3)","e3(1)","e3(2)","e3(3)"'
        line = f.readline()
        while line.strip() != header:
            line = f.readline()
        # Read the connectivity information.
        line = f.readline()
        while line.strip():
            # Read the information, which is comma-delimited.
            num, x1, y1, z1, x2, y2, z2, x3, y3, z3 = line.strip().split(",")
            num, x3, y3, z3 = int(num), float(x3), float(y3), float(z3)
            # Record the information.
            xs[num] = x3
            ys[num] = y3
            zs[num] = z3
            line = f.readline()
    assert sorted(xs) == list(range(min(xs), max(xs) + 1))
    directions = pd.DataFrame.from_dict({"x": xs, "y": ys, "z": zs})
    return directions


def switch_strand(base_nums: Union[Set[int], List[int]], g_ax):
    """
    Return the base numbers corresponding to the other strand.
    :param base_nums: numbers of the bases on the initial strand
    :param g_ax: map of base number to number of paired base
    :return:
    """
    if isinstance(base_nums, set):
        comp_nums = {g_ax[base_num] for base_num in base_nums}
    elif isinstance(base_nums, list):
        comp_nums = [g_ax[base_num] for base_num in base_nums]
    else:
        raise TypeError(type(base_nums))
    return comp_nums


def walk_around_the_block(base_num: int,
                          g_up: Dict[int, Optional[int]],
                          g_dn: Dict[int, Optional[int]],
                          g_ax: Dict[int, Optional[int]],
                          ) -> Tuple[int, int, int, int]:
    """
    Find the number of the base reached by traversing each of the four possible squares.

    :param g_up: map base number to number of the upstream base (or None)
    :param g_dn: map base number to number of the downstream base (or None)
    :param g_ax: map base number to number of the paired base (or None)
    :param base_num: the number of the base to classify
    :return uxux: the base number reached by going upstream, across, upstream, across
    :return xuxu: the base number reached by going across, upstream, across, upstream
    :return dxdx: the base number reached by going downstream, across, downstream, across
    :return xdxd: the base number reached by going across, downstream, across, downstream
    """
    uxux = g_ax.get(g_up.get(g_ax.get(g_up[base_num])))
    xuxu = g_up.get(g_ax.get(g_up.get(g_ax[base_num])))
    dxdx = g_ax.get(g_dn.get(g_ax.get(g_dn[base_num])))
    xdxd = g_dn.get(g_ax.get(g_dn.get(g_ax[base_num])))
    return uxux, xuxu, dxdx, xdxd


def annotate_base(base_num: int,
                  g_up: Dict[int, Optional[int]],
                  g_dn: Dict[int, Optional[int]],
                  g_ax: Dict[int, Optional[int]],
                  ) -> Tuple[str, str, int, bool]:
    """
    Annotate a base in a CanDo file by its strand, structural feature, direction,
    and whether it is opposite or on the same strand as the feature.
    
    :param g_up: map base number to number of the upstream base (or None)
    :param g_dn: map base number to number of the downstream base (or None)
    :param g_ax: map base number to number of the paired base (or None)
    :param base_num: the number of the base to classify
    :return strand: whether base is on SCAF or STAP strand
    :return feature: XO (base participates in crossover or is paired to such a base
                             note that the strand variable does NOT indicate whether
                             the crossover is a scaffold or staple crossover)
                     TM (5' or 3' terminus of a strand, or paired to such a base)
                     EDGE_TM (at the 5' or 3' end of an edge, adjacent to a vertex)
                     VERTEX (an unpaired base in a vertex; currently only applies to staple bases)
                     MIDDLE (none of the above)
    :return direction: for crossover and edge-end, 5 (3) if on 5' (3') side of feature
                       for terminus, 5 (3) if at 5' (3') end of strand
                       for vertex and middle, 0
    :return opposite: True if the feature described is on the opposite strand, else False
    """
    assert base_num is not None
    comp_num = g_ax[base_num]
    assert comp_num != base_num
    strand = STAP if comp_num is None or base_num > comp_num else SCAF
    uxux, xuxu, dxdx, xdxd = walk_around_the_block(base_num, g_up, g_dn, g_ax)
    if base_num == uxux == xuxu == dxdx == xdxd:
        feature, direction, opposite = MIDDLE, 0, False
    elif uxux == xdxd == base_num and xuxu == dxdx and xuxu is not None:
        feature, direction, opposite = XO_2, 5, False
    elif xuxu == dxdx == base_num and uxux == xdxd and uxux is not None:
        feature, direction, opposite = XO_2, 3, False
    elif comp_num is None:  # x
        strand, feature, direction, opposite = STAP, IN_VERTEX, 0, False
    elif g_up[base_num] is None:  # u
        if xdxd is None:
            # 5' terminus is at nick
            feature, direction, opposite = TM5, 0, False
        else: 
            # 5' terminus is adjacent to single crossover
            feature, direction, opposite = TM5_XO, 0, False
    elif g_dn[base_num] is None:  # d
        if xuxu is None:
            # 3' terminus is at nick
            feature, direction, opposite = TM3, 0, False
        else:
            # 3' terminus is adjacent to single crossover
            feature, direction, opposite = TM3_XO, 0, False
    elif g_up[comp_num] is None:  # xu
        if g_dn[g_ax[g_dn[base_num]]] is None: # dxd
            # complementary to 5' terminus at nick
            feature, direction, opposite = TM5, 0, True
        else:
            # complementary to 5' terminus adjacent to single crossover
            feature, direction, opposite = TM5_XO, 0, True
    elif g_dn[comp_num] is None:  # xd
        if g_up[g_ax[g_up[base_num]]] is None: # uxu
            # complementary to 3' terminus is at nick
            feature, direction, opposite = TM3, 0, True
        else:
            # complementary to 3' terminus adjacent to single crossover
            feature, direction, opposite = TM3_XO, 0, True
    elif g_ax[g_up[base_num]] is None:  # ux
        assert strand == STAP
        feature, direction, opposite = VERTEX, 3, False
    elif g_ax[g_dn[base_num]] is None:  # dx
        assert strand == STAP
        feature, direction, opposite = VERTEX, 5, False
    elif g_ax[g_up[comp_num]] is None:  # xux
        assert strand == SCAF
        feature, direction, opposite = VERTEX, 5, False
    elif g_ax[g_dn[comp_num]] is None:  # xdx
        assert strand == SCAF
        feature, direction, opposite = VERTEX, 3, False
    elif g_up[g_ax[g_up[base_num]]] is None:  # uxu
        # the base lies diagonal to a 5' terminus
        # thus its partner must be a 3' terminus or 5' crossover
        # its partner cannot be a 3' terminus b/c then g_dn[comp_num] is None
        # thus its partner must be a 5' crossover
        feature, direction, opposite = XO_1, 5, True
    elif g_dn[g_ax[g_dn[base_num]]] is None:  # dxd
        # the base lies diagonal to a 3' terminus
        # thus its partner must be a 5' terminus or 3' crossover
        # its partner cannot be a 5' terminus or else g_dn[comp_num] is None
        # thus its partner must be a 3' crossover
        feature, direction, opposite = XO_1, 3, True
    elif xuxu is None:  # xuxu
        # the base lies immediately 5' of a 5' terminus
        # thus the base must be a 3' terminus or a 5' crossover
        # it cannot be a 3' terminus or else g_dn[base_num] is None
        # thus the base must be a 5' crossover
        feature, direction, opposite = XO_1, 5, False
    elif xdxd is None:  # xdxd
        # the base lies immediately 3' of a 3' terminus
        # thus the base must be a 5' terminus or a 3' crossover
        # it cannot be a 5' terminus or else g_up[base_num] is None
        # thus the base must be a 3' crossover
        feature, direction, opposite = XO_1, 3, False
    elif uxux is None:  # uxux
        # g_up[g_ax[g_up[base_num]]] must be a vertex base
        # thus g_up[base_num] must be the 3' end of an edge on the scaffold strand
        # thus the base must be the 5' end of another edge on the scaffold strand
        # thus it should be that g_dn[comp_num] is a vertex base and g_ax[g_dn[comp_num]] is None
        # thus this if statement should never be True
        assert False
    elif dxdx is None:  # dxdx
        # g_dn[g_ax[g_dn[base_num]]] must be a vertex base
        # thus g_dn[base_num] must be the 5' end of an edge on the scaffold strand
        # thus the base must be the 3' end of another edge on the scaffold strand
        # thus it should be that g_up[comp_num] is a vertex base and g_ax[g_up[comp_num]] is None
        # thus this if statement should never be True
        assert False
    else:
        # the base should have been annotated by now
        assert False
    assert feature is not None and direction is not None and opposite is not None
    return strand, feature, direction, opposite


def annotate_bases(g_up: Dict[int, Optional[int]],
                   g_dn: Dict[int, Optional[int]],
                   g_ax: Dict[int, Optional[int]],
                   ) -> Dict[Tuple[str, str, int, bool], int]:
    """
    Annotate all of the bases in a CanDo file.
    See annotate_base for more information.

    :param g_up: map base number to number of the upstream base (or None)
    :param g_dn: map base number to number of the downstream base (or None)
    :param g_ax: map base number to number of the paired base (or None)
    :param base_nums: the numbers of the bases to classify
    :return base_annotations: the annotation for each base
    """
    # List the numbers of all of the bases.
    base_nums = sorted(g_up)
    assert base_nums == sorted(g_dn) == sorted(g_ax)
    base_annotations = defaultdict(set)
    for base_num in base_nums:
        # Annotate the base.
        strand, feature, direction, opposite = annotate_base(base_num, g_up, g_dn, g_ax)
        if opposite:
            # Find the number of the base to which the query base is paired.
            comp = g_ax[base_num]
            assert comp is not None
            if comp in base_annotations:
                # If the base is paired to a base that has also been annotated:
                # Ensure that one is scaffold and the other staple.
                assert sorted([strand, base_annotations[comp]["strand"]]) == [SCAF, STAP]
                # Ensure that the locations and directions match.
                assert feature == base_annotations[comp]["location"]
                assert direction == base_annotations[comp]["direction"]
        # Add the base to the annotations.
        base_annotations[strand, feature, direction, opposite].add(base_num)
    # Convert to dict.
    return dict(base_annotations)


'''
def get_base_nums_by_annotations(base_annotations, strands=None, features=None, directions=None, opposites=None):
    """

    :param base_annotations:
    :param strands:
    :param features:
    :param directions:
    :param opposites:
    :return:
    """
    # Assign defaults to missing values.
    strands = [SCAF, STAP] if strands is None else strands
    features = [XO_1, XO_2, EDGE_TM, MIDDLE, TM5, TM3,
                TM_XO, VERTEX] if features is None else features
    directions = [5, 3] if directions is None else directions
    opposites = [True, False] if opposites is None else opposites
    # Retrieve base numbers matching the annotations.
    base_nums = {base_num for annotation in itertools.product(strands, features, directions, opposites) for base_num in base_annotations[annotation]}
    return base_nums
'''


def get_staples_bases_nums(base_annotations_groups: Dict[int, List[int]], g_dn):
    """ Get the residue numbers (according to CanDo numbering) in all of the staples """
    # Find all of the 5' termini of staples.
    termini_5p = set(base_annotations_groups[STAP, STAP_TM5, 0] |
                base_annotations_groups[STAP, STAP_TM5_XO, 0])
    staples_bases_nums = list()
    for terminus_5p in termini_5p:
        # Each terminus corresponds to one staple.
        staple_base_nums = [terminus_5p]
        is_terminus_3p = False
        # Advance through all the base numbers in the staple until reaching the 3' terminus.
        while not is_terminus_3p:
            next_base = g_dn[staple_base_nums[-1]]
            if next_base:
                staple_base_nums.append(next_base)
            else:
                is_terminus_3p = True
        # Add the staple's base numbers to the collection of staples.
        staples_bases_nums.append(staple_base_nums)
    return staples_bases_nums


def get_staples_seqs(staples_bases_nums: List[List[int]], g_ax, scaffold_seq: str) -> List[str]:
    """
    Infer the sequences of staples in a CanDo file based on which base of the scaffold they are paired to
    :param staples_bases_nums: list of staples, each staple represented as a list of its base numbers in CanDo numbering
    :param g_ax: map each base number to the number of the base it is paired to (or None)
    :param scaffold_seq: sequence of the scaffold
    :return: staple_seqs: list of sequences of the staples in the same order as staples_bases
    """
    unpaired = "T"  # unpaired staple residues are T
    staple_seqs = ["".join([seq_utils.comp_base_dna[scaffold_seq[g_ax[base_num] - 1]]
                            if g_ax[base_num] else unpaired for base_num in bases_nums])
                   for bases_nums in staples_bases_nums]
    return staple_seqs
