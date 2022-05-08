import argparse
from collections import defaultdict
import logging

from django.conf import settings
from django.db import connection

logger = logging.getLogger('default')
logger.setLevel(logging.DEBUG)

PLUS = 5
MINUS = 5
REVERSE = True
LENGTH = 20
SEQ_TO_SEARCH = 'GG'
SEQ_REV_TO_SEARCH = 'CC'
MOTIF_LENGTH = 3
# difference between motif length and sequence to search length
# eg. NGG, GG => 1
MOTIF_LENGTH_MINUS = 1
ACGT_REVERSE = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

CAS9 = 1
DCAS9 = 2


def search_sequence(
        full_sequence: str, short_sequence: str, position: int = 0) -> list:
    positions = []
    while True:
        try:
            new_position = full_sequence[position:].index(
                short_sequence) + position
            positions.append(new_position)
            position = new_position + len(short_sequence)
        except ValueError:
            return positions


def search_sequence2(
        full_sequence: str, short_sequence: str, position: int,
        max_position: int) -> list:
    """ Recurently search for sequence with step = 1"""
    #print('full_sequence, short_sequence, position, max_position')
    positions = []
    new_position = None
    while True:
        try:
            new_position = full_sequence[position:].index(short_sequence) \
                           + position
        except ValueError:
            return positions
        if new_position > max_position:
            return positions
        positions.append(new_position)
        position = new_position + 1


def create_sequences(
        full_sequence, length, position, extended=None, mode='cas9'):
    """ Create all sequences that not contain CCN or NGG motif, but 20
    nucleotides next to it (NGG - before NGG found, CCN - after CCN found
    NGG and CGG must be max 5 nucleotides after/before found sequence
    and may overlap with the found sequence"""
    sequences = []
    cas9_from_gg_plus = PLUS
    cas9_to_gg_plus = PLUS-(MOTIF_LENGTH-MOTIF_LENGTH_MINUS)
    dcas9_from_gg_plus = +2
    dcas9_to_gg_plus = LENGTH#+1

    cas9_from_cc_plus = -MINUS
    cas9_to_cc_plus = -MINUS-MOTIF_LENGTH+MOTIF_LENGTH_MINUS
    dcas9_from_cc_plus = -LENGTH - 2
    dcas9_to_cc_plus = -4 #-1

    e_plus = 30 if extended else 0 # FOR EXTENDED

    gg_positions = search_sequence2(
        extended if extended else full_sequence,
        SEQ_TO_SEARCH,
        position+e_plus+(
            cas9_from_gg_plus if mode == 'cas9' else dcas9_from_gg_plus),
        max_position=position+length+e_plus+(
            cas9_to_gg_plus if mode == 'cas9' else dcas9_to_gg_plus))
    for gg_position in set(gg_positions):
        from_seq = gg_position - LENGTH - MOTIF_LENGTH_MINUS-e_plus
        to_seq = gg_position - MOTIF_LENGTH_MINUS-e_plus
        sequence_plus_to_score = None
        result_sequence = None
        if not extended:
            if from_seq < 0 or to_seq > len(full_sequence):
                continue
            else:
                result_sequence = full_sequence[from_seq:to_seq]
            if from_seq-4 < 0 or to_seq+6 > len(full_sequence):
                sequence_plus_to_score = ''
            else:
                sequence_plus_to_score = full_sequence[from_seq-4:to_seq+6]
        else:
            sequence_plus_to_score = extended[from_seq-4+30:to_seq+30+6]
            result_sequence = extended[from_seq+30:to_seq+30]
        sequences.append((result_sequence, '+',
                          gg_position-1-e_plus, sequence_plus_to_score,
                          from_seq, to_seq))

    cc_positions = search_sequence2(
        extended if extended else full_sequence,
        SEQ_REV_TO_SEARCH,
        max(0, position+e_plus+(
            cas9_from_cc_plus if mode == 'cas9' else dcas9_from_cc_plus)),
        max_position=position+length+e_plus+(
            cas9_to_cc_plus if mode == 'cas9' else dcas9_to_cc_plus))
    for cc_position in set(cc_positions):
        from_seq = cc_position + MOTIF_LENGTH-e_plus
        to_seq = cc_position + LENGTH + MOTIF_LENGTH-e_plus
        sequence_plus_to_score = None
        result_sequence = None
        if not extended:
            if from_seq < 0 or to_seq > len(full_sequence):
                continue
            else:
                result_sequence = full_sequence[from_seq:to_seq]
            if from_seq-6 < 0 or to_seq+4 > len(full_sequence):
                sequence_plus_to_score = ''
            else:
                sequence_plus_to_score = full_sequence[from_seq-6:to_seq+4]
        else:
            sequence_plus_to_score = extended[from_seq-6+30:to_seq+4+30]
            result_sequence = extended[from_seq+30:to_seq+30]
        sequence_plus_to_score = ''.join(
            [ACGT_REVERSE[x] for x in reversed(sequence_plus_to_score)])
        # THIS MAY CHANGE IF GUIDE IS CHANGED
        sequences.append(
            (reverse_complementary(result_sequence),
             '-', cc_position-e_plus, sequence_plus_to_score,
             from_seq, to_seq))
    return sequences


def reverse_complementary(sequence):
    """ Create reverse-complementary sequence"""
    new_sequence = ''.join([x for x in reversed(sequence)])
    new_sequence = new_sequence.replace('A', 'Z').replace('T', 'A') \
        .replace('Z', 'T').replace('C', 'Z').replace('G', 'C').replace('Z', 'G')
    return new_sequence


def find_fragments(sequence, motifs, binding_site,
                   binding_site_start, binding_site_end, chromosome,
                   extended, variant):
    sequence = sequence.upper().strip()
    guides_to_create_data = []
    motif_data = {}
    for motif in motifs:
        motif = motif.upper().strip()
        positions = search_sequence(sequence, motif)
        positions.extend(search_sequence(
            sequence, reverse_complementary(motif)))
        for position in sorted(list(set(positions))):
            position_on_chromosome = \
                int(binding_site_start) + position \
                    if chromosome else None

            end_position = position + len(motif)
            motif_data_key = \
                f"{position + 1}_{sequence[position:end_position]}_" \
                f"{position_on_chromosome}_{chromosome}_{binding_site}_" \
                f"{binding_site_start if binding_site_start else None}_" \
                f"{binding_site_end if binding_site_end else None}"
            if not motif_data_key in motif_data:
                motif_data[motif_data_key] = {
                    "query_sequence": motif,
                    "position": position + 1,
                    "found_sequence": sequence[position:end_position],
                    "position_on_chromosome": position_on_chromosome,
                    "chromosome": chromosome,
                    "binding_site_name": binding_site,
                    "binding_site_start": binding_site_start
                    if binding_site_start else None,
                    "binding_site_end": binding_site_end
                    if binding_site_end else None}
            else:
                motif_data[motif_data_key]["query_sequence"] \
                    += f", {motif}"

            guides = create_sequences(
                sequence, len(motif), position, extended,
                ('cas9' if variant == CAS9 else 'dcas9'))

            for guide in guides:
                guides_to_create_data.append(
                    {"found_motif": motif_data_key,  # TODO - change key
                     "position": guide[2] + 1,
                     "sequence": guide[0],
                     "forward_backward": guide[1],
                     "chromosome": chromosome,
                     "sequence_plus_to_score": guide[3],
                     "start_position": guide[4],
                     "end_position": guide[5]})
    return motif_data, guides_to_create_data


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-sequence", help="Sequence file", required=True)
    parser.add_argument("-motifs", help="Motifs file", required=True)
    parser.add_argument("-binding_site", help="Name of the binding site", required=True)
    parser.add_argument("-binding_start", help="Binding site start", required=True)
    parser.add_argument("-binding_end", help="Binding site end", required=True)
    parser.add_argument("-variant", help="Variant (1 for cas9, 2 for dcas9)", required=True)
    parser.add_argument("-extended", help="Extended sequence file", required=False)
    parser.add_argument("-chromosome", help="Chromosome", required=False)
    args = parser.parse_args()

    sequence = None
    extended = None
    motifs = []
    with open(args.sequence, 'r', encoding='utf-8') as f:
        sequence = ''.join(f.readlines()).replace('\n', '')
    if args.extended:
        with open(args.extended, 'r', encoding='utf-8') as f:
            extended = ''.join(f.readlines()).replace('\n', '')
    with open(args.motifs, 'r', encoding='utf-8') as f:
        motifs = f.read().replace('\n', ',').replace(',,', ',').split(',')
    if motifs[-1] == '':
        motifs.pop()

    result = find_fragments(sequence, motifs, args.binding_site, args.binding_start, 
                   args.binding_end, args.chromosome, extended, args.variant)
    print('Motif data', result[0], '\n')
    print('Guides data', result[1])

    # example python3.8 find_motif_in_sequence.py -variant 1 -chromosome chr2 -sequence sequence.txt -motifs motifs.txt -extended extended.txt -binding_site site1 -binding_start 264308 -binding_end 264501




