from Bio.Seq import Seq
from Bio.motifs.matrix import FrequencyPositionMatrix

degenerate_nucleotide = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "AC": "M",
    "AG": "R",
    "AT": "W",
    "CG": "S",
    "CT": "Y",
    "GT": "K",
    "ACG": "V",
    "ACT": "H",
    "AGT": "D",
    "CGT": "B",
    "ACGT": "N",
}
reverse_degenerate_nucleotide = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "N": "N",
}


class CustomFrequencyPositionMatrix(FrequencyPositionMatrix):

    def copy_from_frequency_position_matrix(self, fpm):
        self.length = fpm.length
        self.alphabet = fpm.alphabet
        for letter in fpm.alphabet:
            self[letter] = fpm[letter]

    def degenerate_consensus_cutoff(self, cutoff):

        def get(nucleotide):
            return self[nucleotide][i]

        cutoff = cutoff / 100
        sequence = ""
        for i in range(self.length):
            nucleotides = sorted(self, key=get, reverse=True)
            counts = [self[c][i] for c in nucleotides]
            if counts[3] > cutoff:
                key = "ACGT"
            elif counts[2] > cutoff:
                key = "".join(sorted(nucleotides[:3]))
            elif counts[1] > cutoff:
                key = "".join(sorted(nucleotides[:2]))
            else:
                key = nucleotides[0]
            nucleotide = degenerate_nucleotide.get(key, key)
            sequence += nucleotide
        return Seq(sequence)

    def degenerate_consensus_sum_gte(self, gte_value):

        def get(nucleotide):
            return self[nucleotide][i]

        gte_value = gte_value / 100

        sequence = ""
        for i in range(self.length):
            nucleotides = sorted(self, key=get, reverse=True)
            counts = [self[c][i] for c in nucleotides]
            if counts[0] >= gte_value:
                key = nucleotides[0]
            elif sum(counts[:2]) >= gte_value:
                key = "".join(sorted(nucleotides[:2]))
            elif sum(counts[:3]) >= gte_value:
                key = "".join(sorted(nucleotides[:3]))
            else:
                key = "ACGT"
            nucleotide = degenerate_nucleotide.get(key, key)
            sequence += nucleotide
        return Seq(sequence)

    @staticmethod
    def convert_to_reverse_consensus(sequence):
        return "".join(
            [reverse_degenerate_nucleotide[x] for x in sequence][::-1])
