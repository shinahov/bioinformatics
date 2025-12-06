import collections
from collections import defaultdict

valid_nucleotides = {'A', 'C', 'G', 'T'}
pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

codon_table = {
    "UUU": "F", "UUC": "F",
    "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y",
    "UAA": "Stop", "UAG": "Stop",
    "UGU": "C", "UGC": "C",
    "UGA": "Stop",
    "UGG": "W",

    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H",
    "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",

    "AUU": "I", "AUC": "I", "AUA": "I",
    "AUG": "M",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N",
    "AAA": "K", "AAG": "K",
    "AGU": "S", "AGC": "S",
    "AGA": "R", "AGG": "R",

    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
}


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

def prob_of_dominant_phenotype(p3, p2, p1):
    total_offspring = p1 + p2 + p3
    aaxaa = (p1 / total_offspring) * (p1 - 1) / (total_offspring - 1)
    Aaxaa = (p1 / total_offspring) * (p2 / (total_offspring - 1))
    Aaxaa += (p2 / total_offspring) * (p1 / (total_offspring - 1))
    AaxAa = (p2 / total_offspring) * (p2 - 1) / (total_offspring - 1)
    prob = 1 - (aaxaa + 0.5 * Aaxaa + 0.25 * AaxAa)
    return prob

def translate_rna(rna_sequence):
    protein_sequence = []
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        acid = codon_table.get(codon)
        if acid == "Stop" or acid is None:
            continue
        protein_sequence.append(acid)
    return ''.join(protein_sequence)

def substring_indices(main_string, sub_string):
    indices = []
    for i in range(len(main_string)- len(sub_string)):
        #print(main_string[i:len(sub_string)+i])
        if(main_string[i:len(sub_string)+i] == sub_string):
            indices.append(i+1)
    return ''.join(str(i) + " " for i in indices)

def most_common_ancestor(data):
    data = data.split('>')[1:]
    res = []
    for line in data:
        lines = line.split('\n')
        name = lines[0]
        seq = ''.join(lines[1:])
        res.append((name, seq))
    dna = []
    seq = []
    A = []
    C = []
    G = []
    T = []
    consensus = ""
    for line in res:
        dna.append(line[1])
    for i in range(len(dna[0])):
        for j in range(len(dna)):
            seq.append(dna[j][i])
        Acount = seq.count('A')
        Ccount = seq.count('C')
        Gcount = seq.count('G')
        Tcount = seq.count('T')
        A.append(Acount)
        C.append(Ccount)
        G.append(Gcount)
        T.append(Tcount)
        consensus += max(('A', Acount), ('C', Ccount), ('G', Gcount), ('T', Tcount), key=lambda x: x[1])[0]
        seq = []
        Acount = 0
        Ccount = 0
        Gcount = 0
        Tcount = 0
    A = "A: " + ' '.join(str(i) for i in A)
    C = "C: " + ' '.join(str(i) for i in C)
    G = "G: " + ' '.join(str(i) for i in G)
    T = "T: " + ' '.join(str(i) for i in T)
    return consensus + "\n" + A + "\n" + C + "\n" + G + "\n" + T

def mortal_rabbits(months, months_to_live):
    rabbits = [0] * months_to_live
    rabbits[0] = 1
    for month in range(1, months):
        offsprings = sum(rabbits[1:])
        rabbits = [0] + rabbits[:-1]
        rabbits[0] += offsprings
    return sum(rabbits)

def find_overlaps(data):
    seq = {}
    eqs = defaultdict(list)
    eqp = defaultdict(list)
    name = ""
    for lines in data.splitlines():
        if lines.startswith('>'):
            name = lines[1:]
            seq[name] = ""
        else:
            seq[name] += lines.strip()
    for s in seq:
        suf = seq[s][-3:]
        pref = seq[s][:3]
        eqs[suf].append(s)
        eqp[pref].append(s)
    overlaps = []
    for suf in eqs:
        if suf in eqp:
            for a in eqs[suf]:
                for b in eqp[suf]:
                    if a != b:
                        overlaps.append((a, b))
    return overlaps










