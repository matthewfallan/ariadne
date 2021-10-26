"""
ARIADNE - Vis module

Functions for processing the ASCII visualization of origami structure.
"""

from collections import defaultdict, deque
import re

import seq_utils
import terms


def annotate_bases(vis_file, version):
    """
    Find base annotations from the ASCII visualization file.
    :return:
    """
    with open(vis_file) as f:
        contents = f.read()
    blocks = re.split("Vertex [0-9]+ \(top\) to [0-9]+ \(bottom\), Edge [0-9]+", contents)[1:]
    indexes = list()
    base_annotations = defaultdict(set)
    edges = list()
    for block in blocks:
        double_helices = [deque(), deque()]  # two double helices per edge
        block_lines = [line for line in block.strip().split("\n")[1: -1] if line.strip()]
        n_lines = len(block_lines)
        for i_line, line in enumerate(block_lines, start=1):
            assert len(line) == 26
            idx1 = int(line[0: 6])
            base1 = line[7]
            topo1 = line[9: 13]
            topo2 = line[13: 17]
            base2 = line[18]
            idx2 = int(line[20: 26])
            assert base1 in seq_utils.BASES and base2 in seq_utils.BASES
            indexes.append(idx1)
            indexes.append(idx2)
            num1 = idx1 + 1
            num2 = idx2 + 1
            double_helices[0].appendleft(num1)  # add to 5' side
            double_helices[1].append(num2)  # add to 3' side
            annot1, annot2 = None, None
            if i_line == 1:
                annot1 = terms.EDGE_TM, 3  # at the 3' end of an edge
                annot2 = terms.EDGE_TM, 5  # at the 5' end of an edge
            elif i_line == n_lines:
                annot1 = terms.EDGE_TM, 5  # at the 5' end of an edge
                annot2 = terms.EDGE_TM, 3  # at the 3' end of an edge
            else:
                if version == 1:
                    if topo1 == "| +>":
                        assert topo2 == ">+ |"
                        if block_lines[i_line][9:13] == "| . ": 
                            # next nucleotide is staple end, not staple crossover
                            # block_lines[i_line] is the next line, because i_line = 1 for block_lines[0] 
                            annot1 = terms.STAP_XO_1, 5  # lies just 5' of staple crossover (single; mesojeunction)
                            annot2 = terms.STAP_XO_1, 3  # lies just 3' of staple crossover (single; mesojunction)
                        elif block_lines[i_line][9:13] == "| +<": # next nucleotide is an adjacent staple crossover
                            # next nucleotide is staple end, not staple crossover
                            # block_lines[i_line] is the next line, because i_line = 1 for block_lines[0]    
                            annot1 = terms.STAP_XO_2, 5  # lies just 5' of staple crossover (double; classic junction)
                            annot2 = terms.STAP_XO_2, 3  # lies just 3' of staple crossover (double; classic junction)
                    elif topo2 == ">+ |":
                        raise ValueError()
                    elif topo1 == "| +<":
                        assert topo2 == "<+ |"
                        if block_lines[i_line-2][9:13] == "| v ": # previous nucleotide is a staple end, not staple xover
                            # next nucleotide is staple end, not staple crossover
                            # block_lines[i_line-2] is the previous line, because i_line = 1 for block_lines[0]                            annot1 = terms.STAP_XO_1, 3  # lies just 3' of staple crossover (single; mesojeunction)
                            annot1 = terms.STAP_XO_1, 3  # lies just 3' of staple crossover (single; mesojunction)
                            annot2 = terms.STAP_XO_1, 5  # lies just 5' of staple crossover (single; mesojeunction)
                        elif block_lines[i_line-2][9:13] == "| +>": # previous nucleotide is an adjacent staple crossover
                            # next nucleotide is staple end, not staple crossover
                            # block_lines[i_line-2] is the previous line, because i_line = 1 for block_lines[0]                            annot1 = terms.STAP_XO_2, 3  # lies just 3' of staple crossover (double; classic junction)
                            annot1 = terms.STAP_XO_2, 3  # lies just 3' of staple crossover (double; classic junction)
                            annot2 = terms.STAP_XO_2, 5  # lies just 5' of staple crossover (double; classic junction)
                    elif topo2 == "<+ |":
                        raise ValueError()
                elif version == 2:
                    if topo1 == ">>|>":
                        annot1 = terms.SCAF_XO, 5  # lies just 5' of scaffold crossover
                    if topo2 == ">|>>":
                        annot2 = terms.SCAF_XO, 3  # lies just 3' of scaffold crossover
                    if topo1 == "<<|<":
                        annot1 = terms.SCAF_XO, 3  # lies just 3' of scaffold crossover
                    if topo2 == "<|<<":
                        annot2 = terms.SCAF_XO, 5  # lies just 5' of scaffold crossover
                else:
                    raise ValueError(version)
                if not (annot1 and annot2):
                    if topo1 == "| | ":
                        annot1 = terms.MIDDLE  # normal helix
                    elif topo1 == "| . ":
                        annot1 = terms.STAP_TM, 5  # 5' terminus of staple (not adjacent to single xover)
                    elif topo1 == "| v ":
                        annot1 = terms.STAP_TM, 3 # 3' terminus of staple (not adjacent to single xover)
                    elif topo1 == ". | ":
                        annot1 = terms.SCAF_TM, 3  # 5' terminus of scaffold
                    elif topo1 == "^ | ":
                        annot1 = terms.SCAF_TM, 3  # 3' terminus of scaffold
                    if topo2 == " | |":
                        annot2 = terms.MIDDLE  # normal helix
                    elif topo2 == " . |":
                        annot2 = terms.STAP_TM, 5  # 5' terminus of staple (not adjacent to single xover)
                    elif topo2 == " ^ |":
                        annot2 = terms.STAP_TM, 3  # 3' terminus of staple (not adjacent to single xover)
                    elif topo2 == " | .":
                        annot2 = terms.SCAF_TM, 5  # 5' terminus of scaffold
                    elif topo2 == " | v":
                        annot2 = terms.SCAF_TM, 3  # 3' terminus of scaffold
            assert annot1 and annot2
            base_annotations[annot1].add(num1)
            base_annotations[annot2].add(num2)
        edges.append(double_helices)
    assert sorted(indexes) == list(range(min(indexes), max(indexes))) + [max(indexes)]
    all_nums = sorted([num for annot, nums in base_annotations.items() for num in nums])
    assert all_nums == list(range(all_nums[0], all_nums[-1] + 1))
    edges = [sorted([list(helix1), list(helix2)], key=min) for helix1, helix2 in edges]
    edges = sorted([[helix1, list(reversed(helix2))] for helix1, helix2 in edges],
                   key=lambda x: min(map(min, x)))
    # Convert to dict.
    return dict(base_annotations), edges
