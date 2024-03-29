"""
ARIADNE - ChimeraX module

Interface with UCSF ChimeraX.
"""


import terms


def color_annotations_groups_pdb(base_annotations, cando_num_to_pdb_chain_num, pdb_file, script_file_name=None):
    """

    :param base_annotations_groups:
    :return:
    """
    default_color = "#dddddd"  # light gray
    color_map = {
        (terms.SCAF, terms.SCAF_XO): "#008888",  # dark cyan
        (terms.STAP, terms.SCAF_XO): "#00ffff",  # cyan
        (terms.SCAF, terms.STAP_XO): "#008800",  # green
        (terms.STAP, terms.STAP_XO): "#00ff00",  # bright green
        (terms.SCAF, terms.SCAF_TM): "#884400",  # brown
        (terms.STAP, terms.SCAF_TM): "#ff8800",  # orange
        (terms.SCAF, terms.STAP_TM): "#ff0000",  # red
        (terms.STAP, terms.STAP_TM): "#ff8888",  # pink
        (terms.SCAF, terms.EDGE_TM): "#0000ff",  # blue
        (terms.STAP, terms.EDGE_TM): "#8888ff",  # light blue
        (terms.STAP, terms.VERTEX): "#ffff00",  # yellow
    }
    cmds = [f"open {pdb_file}", f"color {default_color}"]  # default color
    for (strand, feature, direction), base_nums_cando in base_annotations.items():
        if feature == terms.MIDDLE or len(base_nums_cando) == 0:
            continue
        color = color_map[strand, feature]
        chains = [cando_num_to_pdb_chain_num[base_num][0] for base_num in base_nums_cando]
        base_nums_pdb = [cando_num_to_pdb_chain_num[base_num][1] for base_num in base_nums_cando]
        color_items = "".join([f"#1/{chain}:{base_num}" for chain, base_num in zip(chains, base_nums_pdb)])
        cmds.append(f"color {color_items} {color}")
    cmd = "; ".join(cmds)
    if script_file_name:
        #assert not os.path.exists(script_file_name)
        with open(script_file_name, 'w') as f:
            f.write(cmd)
    return cmd