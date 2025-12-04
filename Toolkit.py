import collections

valid_nucleotides = {'A', 'C', 'G', 'T'}
pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def validate_sequence(sequence):
    sequence = sequence.upper()
    for nucleotide in sequence:
        if nucleotide not in valid_nucleotides:
            return False
    return True

def count_nucleotides(dna: str) -> dict:
    counts = {"A": 0, "C": 0, "G": 0, "T": 0}
    for base in dna:
        if base in counts:
            counts[base] += 1
    return counts


def to_rna(sequence):
    sequence = sequence.upper()
    rna_sequence = sequence.replace('T', 'U')
    return rna_sequence

def complement_dna(sequence):
    sequence = sequence.upper()
    complement_sequence = ''.join(pairs[nucleotide] for nucleotide in sequence)
    return complement_sequence

def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length")
    dist = sum(a != b for a, b in (zip(seq1, seq2)))
    return dist