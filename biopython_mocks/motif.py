from Bio.motifs import Motif


class CustomMotif(Motif):
    def __init__(self, alphabet="ACGT", instances=None, counts=None):
        pass

    def create_from_motif(self, motif):
        from .matrix import CustomFrequencyPositionMatrix
        self.name = motif.name
        self.instances = motif.instances
        self.counts = CustomFrequencyPositionMatrix(
            motif.alphabet, {letter: [0, ] for letter in motif.alphabet})
        # here is the most important part:
        self.counts.copy_from_frequency_position_matrix(motif.counts)
        self.length = motif.length
        self.alphabet = motif.alphabet
        self.pseudocounts = motif.pseudocounts
        self.background = motif.background
        self.mask = motif.mask


class IupacMotif(str):
    """This is just to differentiate from string and Motif"""
