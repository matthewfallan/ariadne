"""
ARIADNE - Main module

ARIADNE is a set of tools to analyze and validate
outputs from the program DAEDALUS.

ARIADNE is named after the the Cretan princess who,
in Greek mythology, lent Theseus a thread to help
him navigate the labyrinth built by Daedalus.
"""


from collections import defaultdict
import os
import shutil
import sys
sys.path.append(os.path.dirname(__file__))

import numpy as np
import pandas as pd

import cando
import chimerax
import daedalus
import dssr
import pdb_tools
import plots
import seq_utils


OUTPUTS_DIRECTORY_NAME = "ariadne"


def analyze_design(design_directory, compute_bond_lengths=True, clobber=False):
    design = os.path.basename(design_directory)
    print(design)
    outputs_directory = os.path.join(design_directory, OUTPUTS_DIRECTORY_NAME)
    if os.path.exists(outputs_directory):
        if clobber or seq_utils.yes_no(
                f"Directory {outputs_directory} exists. Overwrite? "):
            shutil.rmtree(outputs_directory)
        else:
            return
    os.mkdir(outputs_directory)
    # Annotate bases with structural data.
    print("\treading PDB file ...")
    pdb_file, cando_file, vis_file = daedalus.get_files(design_directory)
    # process with DSSR
    dssr_info = None
    try:
        #dssr_info = dssr.analyze_pdb(pdb_file)
        pass
    except:
        print("DSSR analysis failed, skipping ...")
    structure = pdb_tools.get_structure(pdb_file)
    model = pdb_tools.get_model(structure)
    print("\treading CanDo file ...")
    base_nums, base_seq, g_up, g_dn, g_ax = cando.get_connectivity(cando_file)
    pair_directions = cando.get_pair_directions(cando_file)
    print("\tannotating bases ...")
    base_annotations, edges, version = daedalus.annotate_bases(g_up, g_dn, g_ax, vis_file)
    cando_num_to_pdb_chain_num = daedalus.map_cando_num_to_pdb_chain_num(
        model, base_annotations, g_dn, g_ax)
    # Compute bond lengths if desired and compile all base information.
    print("\tcompiling all base information ...")
    base_info = daedalus.assemble_base_info(design, model, base_annotations, cando_num_to_pdb_chain_num, pair_directions, g_up, g_dn, g_ax, base_seq, compute_bond_lengths=compute_bond_lengths)
    print("\twriting results files ...")
    base_info_file = os.path.join(outputs_directory, "base_info.tsv")
    base_info.to_csv(base_info_file, sep="\t", index=False)
    # Make a script to color the PDB by annotation.
    chimerax_script = os.path.join(outputs_directory, "chimerax_color.txt")
    chimerax.color_annotations_groups_pdb(base_annotations, cando_num_to_pdb_chain_num, pdb_file, chimerax_script)
    return edges, g_up, g_dn, g_ax, base_info, dssr_info


def analyze_designs(design_directories, compute_bond_lengths=True):
    designs_succeeded = list()
    designs_failed = list()
    designs_base_info = list()
    designs_dssr_info = defaultdict(list)
    designs = list()
    for design_directory_rel in design_directories:
        design_directory_abs = os.path.abspath(design_directory_rel)
        design = os.path.basename(design_directory_abs)
        designs.append(design)
        try:
            edges, g_up, g_dn, g_ax, base_info, dssr_info = analyze_design(design_directory_abs, compute_bond_lengths=compute_bond_lengths, clobber=True)
        except ZeroDivisionError:
            designs_failed.append(design)
        else:
            designs_succeeded.append(design)
            designs_base_info.append(base_info)
            if dssr_info:
                for section, info in dssr_info.items():
                    info["design"] = [design for i in range(info.shape[0])]
                    designs_dssr_info[section].append(info)
    print("Succeeded:", designs_succeeded)
    print("Failed:   ", designs_failed)
    print("plotting DSSR results ...")
    designs_dssr_info = {section: pd.melt(pd.concat(designs),
                                          id_vars=["design"],
                                          value_vars=[field for field, dtype in dssr.DSSR_SECTIONS[section].items() if dtype in [int, float]],
                                          var_name="parameter",
                                          value_name="value",)
                         for section, designs in designs_dssr_info.items()}
    plots.dssr_analysis(designs_dssr_info)
    print("writing bond length file ...")
    designs_base_info = pd.concat(designs_base_info, axis=0, ignore_index=True)
    designs_base_info_file = os.path.join(".", "base_info.tsv")
    designs_base_info.to_csv(designs_base_info_file, sep="\t", index=False)
    if compute_bond_lengths:
        # Melt the table.
        print("reformatting bond lengths ...")
        designs_base_info_long = pd.melt(designs_base_info,
                                         id_vars=daedalus.BOND_ID_VARS,
                                         value_vars=daedalus.BOND_TYPE_VARS,
                                         var_name="bond type",
                                         value_name="bond length (Å)",)
        print("plotting bond lengths ...")
        plots.bond_length_distribution(designs_base_info_long, design_order=designs)


if __name__ == "__main__":
    designs_directories = sys.argv[1:]
    if designs_directories == ["std"]:
        designs_directories = [
            "/Users/mfa/db/LCBB_Matthew_Molly/1_RNAorigami_design/DAEDALUSX_v1_outputs/LibFig_rPB66",
            "/Users/mfa/db/LCBB_Matthew_Molly/1_RNAorigami_design/DAEDALUSX_v1_outputs/LibFig_rOct66",
            "/Users/mfa/db/LCBB_Matthew_Molly/1_RNAorigami_design/DAEDALUSX_v2_outputs/rPB66_v2",
            "/Users/mfa/db/LCBB_Matthew_Molly/1_RNAorigami_design/DAEDALUSX_v2_outputs/rOct66_v2",
        ]
#    analyze_designs(designs_directories)
