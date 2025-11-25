import collections

valid_nucleotides = {'A', 'C', 'G', 'T'}
dictionary_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def validate_sequence(sequence):
    sequence = sequence.upper()
    for nucleotide in sequence:
        if nucleotide not in valid_nucleotides:
            return False
    return True

def count_nucleotides(sequence):
    sequence = sequence.upper()
    for nucleotide in sequence:
        if nucleotide in dictionary_counts:
            dictionary_counts[nucleotide] += 1
    return dictionary_counts

def to_rna(sequence):
    sequence = sequence.upper()
    rna_sequence = sequence.replace('T', 'U')
    return rna_sequence

def complement_dna(sequence):
    sequence = sequence.upper()
    complement_sequence = ''.join(pairs[nucleotide] for nucleotide in sequence)
    return complement_sequence