import pytest

from cas_motif_searcher.biopython_mocks.matrix import CustomFrequencyPositionMatrix


@pytest.mark.parametrize(
    "motif_data, consensus", [
        ({'A': [0.9, 0], 'C': [0.1, 0.5], 'G': [0, 0.5], 'T': [0, 0]}, 'AS'),
        ({'A': [0.51, 0], 'C': [0.24, 0.5], 'G': [0.24, 0.5], 'T': [0.01, 0]},
         'AS'),
        ({'A': [0.51, 0], 'C': [0.26, 0.5], 'G': [0.23, 0.5], 'T': [0.0, 0]},
         'MS'),
        ({'A': [0.49, 0], 'C': [0.23, 0.5], 'G': [0.23, 0.5], 'T': [0.05, 0]},
         'NS'),

    ])
def test_degenerate_consensus(motif_data, consensus):
    """If first most frequent nucleotide > 50 % and next < 2 times smaller 1 nt
    elif two nucleodites > 75% - consensus 2, else N"""
    matrix = CustomFrequencyPositionMatrix(alphabet='ACGT', values=motif_data)
    assert matrix.degenerate_consensus == consensus


@pytest.mark.parametrize(
    "motif_data, cutoff, consensus", [
        ({'A': [0.9, 0], 'C': [0.1, 0.5], 'G': [0, 0.5], 'T': [0, 0]}, 5, 'MS'),
        ({'A': [0.96, 0.25], 'C': [0.04, 0.25], 'G': [0, 0.5], 'T': [0, 0]},
         5, 'AV'),
    ])
def test_degenerate_consensus_cutoff(motif_data, cutoff, consensus):
    matrix = CustomFrequencyPositionMatrix(alphabet='ACGT', values=motif_data)
    assert matrix.degenerate_consensus_cutoff(cutoff) == consensus


@pytest.mark.parametrize(
    "motif_data, gte, consensus", [
        ({'A': [0.9, 0], 'C': [0.1, 0.5], 'G': [0, 0.5], 'T': [0, 0]},
         95, 'MS'),
        ({'A': [0.96, 0.25], 'C': [0.04, 0.25], 'G': [0, 0.5], 'T': [0, 0]},
         95, 'AV'),
        ({'A': [0.94, 0], 'C': [0.02, 0.5], 'G': [0, 0.5], 'T': [0, 0]},
         95, 'MS'),
    ])
def test_degenerate_consensus_sum_gte(motif_data, gte, consensus):
    matrix = CustomFrequencyPositionMatrix(alphabet='ACGT', values=motif_data)
    assert matrix.degenerate_consensus_sum_gte(gte) == consensus


@pytest.mark.parametrize(
    "sequence, reversed", [
        ('CAGAT', 'ATCTG'),
        ('KYSWRMBDHVN', 'NBDHVKYWSRM'), # W i S są torzsame!!!
    ])
def test_convert_to_reverse_consensus(sequence, reversed):
    matrix = CustomFrequencyPositionMatrix(
        alphabet='ACGT',
        values={'A': [0.9, 0], 'C': [0.1, 0.5], 'G': [0, 0.5], 'T': [0, 0]})
    assert matrix.convert_to_reverse_consensus(sequence) == reversed