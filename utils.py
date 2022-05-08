from collections import defaultdict
import linecache
import logging
import os

from django.conf import settings
from django.db import connection

from exceptions import (
    ImproperCoordinateFileException, 
    ImproperNucleotideException, 
    NoSequenceException)


def get_genome_fragment_by_coordinates(genome_type, chromosome, start, end):
    """
    example of proper sequence
    hg19_dna range=chr19:42381139-42381289
GCCACCTCTCAGGGGAATTGTGGCCACCACAGGTGCAGGGAGCAGTTTCT
CTCCACTCACAGCCTGAAGCATACCCGGCAGGGGCTGTCCCCAGGCCCAA
CAAGCAAAGGGCCCAGTAGCGAGGGCCACTGGAGCCCATCTCCGGGGGGC
T
    :param genome_type:
    :param chromosome:
    :param start:
    :param end:
    :return:
    """
    try:
        start = int(start)
        end = int(end)
    except ValueError:
        return (None, None)
    start -= 31
    end += 30
    if start < 0 or end < 0:
        return (None, None)
    #/organisms/hg19_all$ tail -n 300 chr7.fa
    filename = os.path.join(
        os.path.dirname(__file__), 'organisms', f'{genome_type}', f'{chromosome}.fa')
    lines_start = 2 + start // 50
    lines_end = 2 + end // 50
    cutoff_start = start % 50
    cutoff_end = 50 - (end % 50)
    lines = []
    for line_no in range(lines_start, lines_end+1):
        lines.append(linecache.getline(filename, line_no).strip())
    sequence = ''.join(lines)[cutoff_start:-cutoff_end]
    return (sequence[30:-30] if sequence else None), \
           (sequence if sequence else None)


def fix_sequences(sequences, iupac=False):
    """Change all letters in the sequences to uppercase and
    check if all of the sequences contains only ACTG characters"""
    letters = 'ACTG' if not iupac else 'ACGTMRWSYKVHDBN'
    new_sequences = []
    for sequence in sequences:
        if isinstance(sequence, tuple): # for fasta names
            sequence_upper = sequence[0].upper()
            if not set(sequence_upper).issubset(set(letters)):
                raise ImproperNucleotideException(f'Wrong sequence: {sequence[0]}')
            new_sequences.append((sequence_upper.upper(), sequence[1]))
        else:
            sequence_upper = sequence.upper()
            if not set(sequence_upper).issubset(set(letters)):
                raise ImproperNucleotideException(f'Wrong sequence: {sequence}')
            new_sequences.append(sequence_upper.upper())
    return new_sequences


def add_names_to_sequences(
        sequences, chromosome=None, start_position=None, end_position=None):
    # for fasta names
    if len(sequences) > 0 and isinstance(sequences[0], tuple):
        return [(sequence[0], (sequence[1] if sequence[1]
                               else f"{settings.BS_SEQ_NAME}{i+1}"),
                 start_position, end_position, chromosome, None)
                for i, sequence in enumerate(sequences)]

    return [(sequence, (f"{chromosome}_" if chromosome else "")
             + f"{settings.BS_SEQ_NAME}{i+1}",
             start_position, end_position, chromosome, None)
            for i, sequence in enumerate(sequences)]


def process_fasta(lines_generator, iupac=False):
    sequences = []
    current_sequence = []
    name = ''
    for line in lines_generator:
        if line.startswith('>'):
            if current_sequence:
                sequences.append((''.join(current_sequence), name))
                current_sequence = []
            name = line[1:].strip()
        else:
            current_sequence.append(line.strip())
    if current_sequence:
        sequences.append((''.join(current_sequence).upper(), name))
    return fix_sequences(sequences, iupac)


def extract_sequences_from_file(filename, iupac=False):
    with open(filename) as f:
        return process_fasta(f, iupac)


def extract_sequences_from_text_field(content, iupac=False):
    if content.startswith('>'):
        return process_fasta(content.split('\n'), iupac)
    for letter in ",;:|":
        content = content.replace(letter, ' ')
    elements = content.split()
    if isinstance(elements, str):
        elements = [elements, ]
    return fix_sequences(elements, iupac)



def process_coordinates_file(query, filename):
    sequences = []
    chromomes_count = defaultdict(int)
    chromosomes_data = defaultdict(list)
    with open(filename) as f:
        i = 1
        for line in f:
            elements = line.strip().split()
            if len(elements) < 3:
                if len(elements) == 0:
                    continue
                raise ImproperCoordinateFileException(
                    f"Too short line: \"{line.strip()}\", "
                    f"use format like: \"chr1 111 222\" for each line")
            if not elements[0].startswith('chr'):
                raise ImproperCoordinateFileException(
                    f"Line does not start with chr: \"{line.strip()}\", "
                    f"use format like: \"chr1 111 222\" for each line")
            try:
                int(elements[1])
            except ValueError:
                raise ImproperCoordinateFileException(
                    f"Improper first coordinate: \"{line.strip()}\", "
                    f"use format like: \"chr1 111 222\" for each line")
            try:
                int(elements[2])
            except ValueError:
                raise ImproperCoordinateFileException(
                    f"Improper second coordinate: \"{line.strip()}\", "
                    f"use format like: \"chr1 111 222\" for each line")

            chromosomes_data[elements[0]].append(elements)

    for chromosome, all_elements in chromosomes_data.items():
        linecache.clearcache()
        for elements in all_elements:
            chromomes_count[elements[0]] += 1 #for numeration by the chromosome
            sequence, extended = get_genome_fragment_by_coordinates(
                genome_type=query.genome, chromosome=elements[0],
                start=elements[1], end=elements[2])
            if not sequence:
                raise NoSequenceException(
                    f'Sequence cannot be retrieved ({elements[0]} '
                    f'{elements[1]} {elements[2]})')
            if len(elements) > 3:
                sequences.append(
                    (sequence.upper(), elements[3], elements[1], elements[2],
                     elements[0], extended.upper()))
            else:
                sequences.append(
                    (sequence.upper(),
                     f"{elements[0]}_{settings.BS_SEQ_NAME}"
                     f"{chromomes_count[elements[0]]}",
                     elements[1], elements[2],
                     elements[0], extended.upper()))
            i += 1
    linecache.clearcache()
    return sequences

