"""
DeBruijn Sequence Generator

For generating reverse-complement-free nucleic acid strands using deBruijn graphs.
"""

import argparse
from collections import Counter, defaultdict
import sys
import itertools
import random

import matplotlib.pyplot as plt


bases = "ACGT"


def rc(seq):
    rcs = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(reversed([rcs[base] for base in seq]))


def assemble(kmers):
    return kmers[0] + ''.join([kmer[-1] for kmer in kmers[1:]])


def disassemble(seq, k):
    return [seq[i: i + k] for i in range(len(seq) - k + 1)]


def gc_content(seq):
    return (seq.count('C') + seq.count('G')) / len(seq)


def valid_kmer(kmer):
    return kmer != rc(kmer) and not any([forbidden in kmer for forbidden in forbidden_seqs])


def validate_seq(seq, k, forbidden_seqs=None):
    if forbidden_seqs is None:
        forbidden_seqs = list()
    kmers = Counter(disassemble(seq, k))
    # No forbidden sequences.
    if forbidden_seqs:
        forbidden_seq_list = [forbid for forbid in forbidden_seqs if forbid in seq]
        if forbidden_seq_list:
            return ValueError(f"forbidden sequences: {forbidden_seq_list}")
    # No kmer should appear twice.
    repeated_seq_list = [kmer for kmer, count in kmers.items() if count > 1]
    if repeated_seq_list:
        return ValueError(f"repeated sequences: {repeated_seq_list}")
    # No kmer should have its reverse complement appear.
    rc_kmers = [kmer for kmer in kmers if rc(kmer) in kmers]
    if rc_kmers:
        return ValueError(f"kmers with reverse complements: {rc_kmers}")
    return


def get_random_kmer(k, start=None, max_tries=100):
    valid = False
    tries = 0
    while not valid:
        tries += 1
        if tries > max_tries:
            raise ValueError("Failed to generate valid random kmer")
        if start is None:
            kmer = "".join([random.choice(bases) for i in range(k)])
        else:
            assert len(start) <= k
            kmer = start + "".join([random.choice(bases) for i in range(k - len(start))])
        valid = valid_kmer(kmer)
    return kmer


promoter_starts = {
    "T7": "G",
}


def generate_sequence(seq_len_target, k, forbidden_seqs=None, promoter=None, attempts_limit=1000):
    if forbidden_seqs is None:
        forbidden_seqs = list()
    kmers_next = dict()
    kmers_excl = defaultdict(set)
    
    attempted_seqs = list()
    kmer_stack = list()
    length_max = 0
    attempts = 0
    failed = False
    seq_len = 0
    # Build the sequence until it reaches full length.
    while seq_len < seq_len_target:
        if kmer_stack:
            # If there are already kmers in the sequence:
            # Get the last kmer in the sequence.
            last_kmer = kmer_stack[-1]
            # Find all kmers that can come next, excluding those that already appear and whose reverse complements already appear.
            if last_kmer not in kmers_next:
                kmers_next[last_kmer] = [last_kmer[1:] + base for base in bases if valid_kmer(last_kmer[1:] + base)]
            kmer_options = [kmer for kmer in kmers_next[last_kmer] if kmer not in kmers_excl[last_kmer] and kmer not in kmer_stack and rc(kmer) not in kmer_stack]
            if kmer_options:
                # If there are options, choose one.
                kmer_weights = [relative_weights[kmer[-1]] for kmer in kmer_options]
                kmer_next = random.choices(kmer_options, weights=kmer_weights, k=1)[0]
                # Add the kmer to the sequence.
                kmer_stack.append(kmer_next)
                # Prevent the kmer from being chosen again in the case of backtracking.
                kmers_excl[last_kmer].add(kmer_next)
            else:
                # If not, remove last kmer.
                kmer_stack.pop()
                # Also remove any excluded kmers from following that kmer, in case that kmers is reached again in a different context.
                kmers_excl[last_kmer] = set()
        else:
            # If the sequence is empty, initialize with a kmer.
            if promoter is None:
                kmer_start = None
            else:
                try:
                    kmer_start = promoter_starts[promoter]
                except KeyError:
                    raise ValueError(f"unrecognized promoter: '{promoter}'")
            kmer_next = get_random_kmer(k, start=kmer_start)
            kmer_stack.append(kmer_next)
        seq_len = len(kmer_stack) + k - 1
        print(f"current length: {seq_len}", end="\r")
        if seq_len > length_max:
            length_max = seq_len
            attempts = 0
        else:
            attempts += 1
            if attempts > attempts_limit:
                break
    
    print()
    if len(kmer_stack) == seq_len_target - k + 1:
        sequence = assemble(kmer_stack)
        errors = validate_seq(sequence, k, forbidden_seqs=forbidden_seqs)
        if errors:
            raise errors
        print("success")
        return sequence
    else:
        print("failed")


def analyze_content(seq, window=20):
    print(f"GC content: {round(gc_content(seq), 3)}")
    print(", ".join([f"{base}: {round(seq.count(base) / len(seq), 3)}" for base in bases]))
    gc_content_sliding(seq, window=window)


def gc_content_sliding(seq, window=20):
    gc_content_sliding = list(map(gc_content, disassemble(seq, window)))
    plt.plot(gc_content_sliding)
    plt.xlabel("position")
    plt.ylabel("GC content")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("length", type=int)
    parser.add_argument("k", type=int)
    parser.add_argument("--output", nargs="?", default=None)
    parser.add_argument("--promoter", nargs="?", default=None)
    parser.add_argument("--max_consecutive", nargs="?", default=None, type=int)
    parser.add_argument("--no_cpg", action="store_true")
    parser.add_argument("--no_start", action="store_true")
    parser.add_argument("--gc_target", type=float, default=0.5)
    args = parser.parse_args()

    seq_len_target = args.length
    k = args.k
    fout = args.output
    promoter = args.promoter
    max_cons = args.max_consecutive
    gc = args.gc_target

    forbidden_seqs = list()
    if args.no_start:
        forbidden_seqs.append("ATG")
    if args.no_cpg:
        forbidden_seqs.append("CG")
        relative_weights = {"A": (1-gc)/2, "C": gc/2, "G": 1.5*gc/2, "T": (1-gc)/2}
    else:
        relative_weights = {"A": (1-gc)/2, "C": gc/2, "G": gc/2, "T": (1-gc)/2}
    if max_cons is not None:
        forbidden_seqs.extend([base * (max_cons + 1) for base in bases])

    seq = generate_sequence(seq_len_target, k, forbidden_seqs=forbidden_seqs, promoter=promoter)
    if seq is not None:
        if fout:
            with open(fout, "w") as f:
                f.write(seq)
        else:
            print(seq)
            analyze_content(seq)
