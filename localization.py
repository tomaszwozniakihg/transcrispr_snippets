import argparse
from collections import defaultdict
import logging
import os


logger = logging.getLogger('default')
logger.setLevel(logging.DEBUG)

GENE = 1
INTRON = 2
INTER_GENE = 3
UTR = 4
NOT_SET = 20
LOCALIZATION_MAPPER = {
    GENE:'Coding exon',
    INTRON:'Intron',
    INTER_GENE:'Intergenic',
    UTR:'Non-coding exon',
    NOT_SET:'Not set'}


def get_localizations(genome_name, positions):
    logger.debug('get localizations start ' + str(genome_name) + " | "
                + str(positions))
    ordered_positions = defaultdict(list)
    for position in positions:
        ordered_positions[position[0]].append(
            {'start': position[1], 'end': position[2], 'strand': position[3],
             'annotations': [], 'prev': [], 'next': [], 'current_gene': None,
             'previous_gene': None, 'next_gene': None,
             'previous_gene_start': None, 'next_gene_start': None})

    for chromosome in ordered_positions:
        ordered_positions[chromosome] = sorted(
            ordered_positions[chromosome], key=lambda x: x['start'])

    for chromosome in ordered_positions:
        chromosome_filename = f"{genome_name}_{chromosome}.bed"
        chromosome_filepath = os.path.join(
            os.path.dirname(__file__), 'localization_data', 'elements_coordinates',
            chromosome_filename)
        current_position = 0
        current_position_end = 0
        positions_buffer = []
        all_positions_buffer = []
        with open(chromosome_filepath) as f:
            for position in ordered_positions[chromosome]:
                logger.debug(f'position {position}')
                position_start = position['start']
                position_end = position['end'] -1 #FIX for UCSC notation
                line_elements = None
                current_gene = None

                # clear buffer
                new_buffered = []
                for buffered in all_positions_buffer:
                    if int(buffered[2]) >= position_start:
                        new_buffered.append(buffered)
                all_positions_buffer = new_buffered

                positions_buffer = [x for x in all_positions_buffer
                                    if int(x[1]) <= position_end]

                # find starting place
                while current_position_end < position_start:
                    current_line = f.readline().strip()
                    line_elements = current_line.split()
                    if len(line_elements) < 5:
                        break #there is nothing left
                    line_start = int(line_elements[1])
                    line_end = int(line_elements[2])
                    current_gene = line_elements[3].split(".")[0] \
                        if line_elements[4] == 'coding_exons' else None
                    current_position = line_start
                    current_position_end = line_end
                if line_elements:
                    if current_position <= position_end:
                        positions_buffer.append(line_elements)
                    all_positions_buffer.append(line_elements)
                logger.debug(f'post_finding_starting_place {line_elements} '
                            f'{positions_buffer} {all_positions_buffer}')
                # get all that covers given place
                while current_position <= position_end:
                    logger.debug(f'inside while loop, {current_position}-'
                                 f'{current_position_end}:{position_start}'
                                 f'-{position_end}')
                    # write current gene here, otherwise next line read
                    # overwrites this, and it may be outside given range:
                    if current_gene:
                        position['current_gene'] = current_gene

                    if line_elements and len(line_elements) >= 5:
                        if current_position_end > position_start:
                            logger.debug(f'adding_elements {line_elements}')
                            positions_buffer.append(line_elements)
                    current_line = f.readline().strip()
                    logger.debug('start end loop current line '
                                + str(current_line))
                    if not current_line:
                        break
                    line_elements = current_line.split()
                    if len(line_elements) < 5:
                        break
                    line_start = int(line_elements[1])
                    line_end = int(line_elements[2])
                    current_gene = line_elements[3].split(".")[0] \
                        if line_elements[4] == 'coding_exons' else None
                    current_position = line_start
                    current_position_end = line_end
                    logger.debug(f'next_elements: {line_elements}')
                    if line_elements and len(line_elements) >= 5:
                        all_positions_buffer.append(line_elements)
                for buffered in positions_buffer:
                    if len(buffered) >= 5:
                        position['annotations'].append(buffered[4])
                    else:
                        logger.warning('Wrong buffered data (genomic '
                                       'localization) {}'.format(
                            str(positions_buffer)))
                logger.debug(f'all buffered {positions_buffer}')
                logger.debug(f'buffer for longer {all_positions_buffer}')
    return ordered_positions


def get_gene_start_before_after(genome_name, ordered_positions):
    filepath = os.path.join(
        os.path.dirname(__file__), 'localization_data', 'genes_coordinates',
        f'{genome_name}.bed')
    current_gene_name = None
    previous_gene_name = None
    current_start = -1000000
    current_transcription_start = -10000000
    previous_transcription_start = -10000000
    current_chromosome = None
    with open(filepath) as f:
        for chromosome in sorted(
                list(ordered_positions.keys()),
                key=lambda x: (not x[3:].isdigit(), len(x), x)):
            for position in ordered_positions[chromosome]:
                #print(position)
                while current_chromosome != chromosome:
                    current_line = f.readline().strip()
                    if not current_line:
                        break
                    line_elements = current_line.split()
                    if len(line_elements) < 6:
                        break
                    current_gene_name = line_elements[3]
                    current_start = int(line_elements[1]) \
                        if line_elements[5] == '+' else int(line_elements[2])
                    current_transcription_start = int(line_elements[1]) + 1 \
                        if line_elements[5] == '+' else int(line_elements[2])
                    current_chromosome = line_elements[0]

                while position['start'] > current_start:
                    previous_gene_name = current_gene_name
                    previous_transcription_start = current_transcription_start
                    current_line = f.readline().strip()
                    if not current_line:
                        break
                    line_elements = current_line.split()
                    if len(line_elements) < 6:
                        break
                    current_gene_name = line_elements[3]
                    current_transcription_start = int(line_elements[1]) + 1 \
                        if line_elements[5] == '+' else int(line_elements[2])
                    current_start = int(line_elements[1]) \
                        if line_elements[5] == '+' else int(line_elements[2])
                    current_chromosome = line_elements[0]
                position["previous_gene"] = previous_gene_name
                position["next_gene"] = current_gene_name
                position["previous_gene_start"] = previous_transcription_start \
                    if previous_transcription_start > 0 else None
                position["next_gene_start"] = current_transcription_start \
                    if current_transcription_start > 0 else None
    return ordered_positions


def get_localizations_from_query(motifs_filename, output_filename, genome):
    positions = []
    motifs = []
    with open(motifs_filename, 'r', encoding='utf-8') as f:
        for line in f:
            if line.strip():
                line_elements = line.strip().split(',')
                motifs.append(FoundMotif(*line_elements))

    for motif in motifs:
        if not motif.chromosome:
            continue
        #position_on_chromosome
        coordinates_end = motif.position_on_chromosome + \
                          len(motif.query_sequence.split(',')[0])
        forward_backward = "+"
        positions.append(
            (motif.chromosome, motif.position_on_chromosome, coordinates_end,
             forward_backward))

    locations = get_localizations(genome, positions)
    locations = get_gene_start_before_after(genome, locations)
    for motif in motifs:
        if not motif.chromosome:
            continue
        for chromosome, chromosome_locations in locations.items():
            for location in chromosome_locations:
                if (motif.chromosome == chromosome
                        and location['start'] == motif.position_on_chromosome
                        and location['end'] == motif.position_on_chromosome
                            + len(motif.query_sequence.split(',')[0])):

                    if 'coding_exons' in location['annotations']:
                        motif.motif_location = GENE
                    elif 'utr' in location['annotations']:
                        motif.motif_location = UTR
                    elif 'introns' in location['annotations']:
                        motif.motif_location = INTRON
                    else:
                        motif.motif_location = INTER_GENE
                    motif.current_gene = location['current_gene']
                    motif.previous_gene = location['previous_gene']
                    motif.next_gene = location['next_gene']
                    motif.previous_gene_start = location['previous_gene_start']
                    motif.next_gene_start = location['next_gene_start']

    with open(output_filename, 'w', encoding='utf-8') as f:
        f.write("id,chromosome,position_on_chromosome,query_sequence,motif_location,current_gene,"
                "previous_gene,next_gene,previous_gene_start,next_gene_start\n")
        for motif in motifs:
            motif.write(f)



# code that can be used instead of ORM (database connection)
class FoundMotif:
    current_id = 1

    def __init__(self, chromosome, position_on_chromosome, query_sequence):
        self.id = FoundMotif.current_id
        FoundMotif.current_id += 1
        self.chromosome = chromosome
        self.position_on_chromosome = int(position_on_chromosome)
        self.query_sequence = query_sequence

        self.motif_location = NOT_SET
        self.current_gene = None
        self.previous_gene = None
        self.next_gene = None
        self.previous_gene_start = None
        self.next_gene_start = None

    def write(self, handler):
        handler.write(f"{self.id},{self.chromosome},{self.position_on_chromosome},"
                      f"{self.query_sequence},{LOCALIZATION_MAPPER[self.motif_location]},"
                      f"{self.current_gene},{self.previous_gene},{self.next_gene},"
                      f"{self.previous_gene_start},{self.next_gene_start}\n")



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="Input file", required=True)
    parser.add_argument("-o", help="Output file", required=True)
    parser.add_argument("-g", help="Genome", required=True)
    args = parser.parse_args()
    get_localizations_from_query(args.i, args.o, args.g)
    # use motifs.csv if you wan to check functionality (chromosome,position_on_chromosome,query_sequence)
    #python3.8 localization.py -i motifs.csv -o o.txt -g hg19






