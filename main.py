import Toolkit


if __name__ == '__main__':
    pairs = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    DNA_sequence = ""

    rna = Toolkit.to_rna(DNA_sequence)
    print(rna)



