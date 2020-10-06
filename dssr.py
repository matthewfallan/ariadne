"""
ARIADNE - DSSR module

Interface with 3DNA-DSSR
"""

import os

import numpy as np
import pandas as pd

import seq_utils

DSSR_CMD = "x3dna-dssr"
MISSING = "----"

DSSR_SECTIONS = {
    "Local base-pair parameters": {
        "bp": str,
        "Shear": float,
        "Stretch": float,
        "Stagger": float,
        "Buckle": float,
        "Propeller": float,
        "Opening": float,
    },
    "Local base-pair step parameters": {
        "bp": str,
        "Shift": float,
        "Slide": float,
        "Rise": float,
        "Tilt": float,
        "Roll": float,
        "Twist": float,
    },
    "Local base-pair helical parameters": {
        "bp": str,
        "X-disp": float,
        "Y-disp": float,
        "h-Rise": float,
        "Incl.": float,
        "Tip": float,
        "h-Twist": float,
    },
    "Simple base-pair parameters based on RC8--YC6 vectors": {
        "bp": str,
        "Shear": float,
        "Stretch": float,
        "Stagger": float,
        "Buckle": float,
        "Propeller": float,
        "Opening": float,
    },
    "Simple step parameters based on consecutive C1'-C1' vectors": {
        "bp": str,
        "Shift": float,
        "Slide": float,
        "Rise": float,
        "Tilt": float,
        "Roll": float,
        "Twist": float,
    },
    "Simple helical parameters based on consecutive C1'-C1' vectors": {
        "bp": str,
        "X-disp": float,
        "Y-disp": float,
        "h-Rise": float,
        "Incl.": float,
        "Tip": float,
        "h-Twist": float,
    },
    "Classification of dinucleotide steps in right-handed nucleic acid structures": {
        "bp": str,
        "Xp": float,
        "Yp": float,
        "Zp": float,
        "XpH": float,
        "YpH": float,
        "ZpH": float,
        "Form": str,
        "ABI": float,
    },
}


def parse_dssr_section(output_file, line):
    section = None
    while section is None:
        sections = [section for section in DSSR_SECTIONS if section in line]
        if sections:
            assert len(sections) == 1
            section = sections[0]
        line = output_file.readline()
        if line == "":
            # end of file
            return None
    fields = DSSR_SECTIONS[section]
    field_names = set(fields)
    while set(line.split()) != field_names:
        line = output_file.readline()
        if line == "":
            # end of file
            return None
    field_names_order = line.split()
    field_types = {field: fields[field] for field in field_names_order}
    section_data = {field: list() for field in field_names if field_types[field] in [int, float]}
    line = output_file.readline()
    while line.strip() != "":
        for field, value in zip(field_names_order, line.split()[1:]):
            if field not in section_data:
                continue
            if value == MISSING:
                section_data[field].append(np.nan)
            else:
                section_data[field].append(field_types[field](value))
        line = output_file.readline()
    section_data = pd.DataFrame.from_dict(section_data)
    return section, section_data


def parse_dssr_output(output_file_name):
    with open(output_file_name) as output_file:
        output_data = dict()
        for i in range(len(DSSR_SECTIONS)):
            result = parse_dssr_section(output_file, output_file.readline())
            if result:
                section, section_data = result
                output_data[section] = section_data
            else:
                assert False
    return output_data


def analyze_pdb(pdb_file, clobber=False):
    """
    Determine the results of DSSR.
    :return:
    """
    base_name, ext = os.path.splitext(pdb_file)
    output_file_name = f"{base_name}_dssr.txt"
    if clobber or not os.path.isfile(output_file_name):
        cmd = f"{DSSR_CMD} --analyze -i={seq_utils.escape_path(pdb_file)} -o={seq_utils.escape_path(output_file_name)}"
        xc = os.system(cmd)
        assert xc == 0
    output_data = parse_dssr_output(output_file_name)
    return output_data

