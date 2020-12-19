"""
ARIADNE - Seq Funcs module

Basic tools for working with sequences.
"""

scaffold_chain_id = "A"

BASES_DNA = {"A", "C", "G", "T"}
BASES_RNA = {"A", "C", "G", "U"}
BASES = BASES_DNA | BASES_RNA

comp_base_dna = {"A": "T", "C": "G", "G": "C", "T": "A", "U": "A"}

def get_rc_dna(seq):
    return "".join(map(lambda x: comp_base_dna[x], reversed(seq)))


comp_base_rna = {"A": "U", "C": "G", "G": "C", "T": "A", "U": "A"}

def get_rc_rna(seq):
    return "".join(map(lambda x: comp_base_rna[x], reversed(seq)))


def switch_direction(direction):
    return {"0": "0", "5": "3", "3": "5", 0: 0, 5: 3, 3: 5}[direction]


def yes_no(question):
    yes = ["y", "ye", "yes"]
    no = ["n", "no"]
    answer = None
    while answer not in yes + no:
        answer = input(question)
    return answer in yes


def escape_path(path):
    escaped = path.replace("(", "\(").replace(")", "\)").replace(" ", "\ ")
    return escaped


def get_gc_content(seq):
    return (seq.count("G") + seq.count("C")) / len(seq)