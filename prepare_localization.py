import argparse
from collections import defaultdict
import os



def create_directories():
    base_dir = os.path.join(os.path.dirname(__file__), 'localization_data')
    dirs_to_be_created = [
        base_dir, os.path.join(base_dir, 'temp'), 
        os.path.join(base_dir, 'genes_coordinates'),
        os.path.join(base_dir, 'elements_coordinates')]
    for directory in dirs_to_be_created:
        if not os.path.exists(directory):
            os.mkdir(directory)

    
def download_single_genome_data(genome_name):
    """This function download data from UCSC and save it to specific file,
    then sort it
    more details on how does it work:
    https://genome.ucsc.edu/cgi-bin/hgTables
    http://genomewiki.ucsc.edu/index.php/Programmatic_access_to_the_Genome_Browser

    "mysql hg38  -h genome-mysql.soe.ucsc.edu -u genome -A -e 'select * from ncbiRefSeqCurated' -NB > out.txt"

    file format:
    0	NM_001276352.2	chr1	-	67092164	67134970	67093579	67127240	9	67092164,67096251,67103237,67111576,67115351,67125751,67127165,67131141,67134929,	67093604,67096321,67103382,67111644,67115464,67125909,67127257,67131227,67134970,	0	C1orf141	cmpl	cmpl	2,1,0,1,2,0,0,-1,-1,
    """
    genome_dir = os.path.join(
            os.path.dirname(__file__), 'localization_data', 'temp')
    genome_file = os.path.join(genome_dir, f'{genome_name}.txt')
    command = f"mysql {genome_name} -h genome-mysql.soe.ucsc.edu -u genome " \
              f"-A -e 'select * from ncbiRefSeqCurated' -NB > {genome_file}"
    os.system(command)
    sorted_genome_file = os.path.join(genome_dir, f'{genome_name}_sorted.txt')

    with open(genome_file, 'r') as f_in:
        with open(sorted_genome_file, 'w') as f_out:
            lines = [line.strip().split() for line in f_in.readlines()]
            lines = [(*l[:4], int(l[4]), *l[5:]) for l in lines]
            lines = sorted(
                lines,
                key=lambda x: (not x[2][3:].isdigit(), len(x[2]), x[2], x[4]))
            content = "\n".join(["\t".join([str(x) for x in line])
                                 for line in lines])
            f_out.write(content)
    os.remove(genome_file)


def create_gene_and_elements_files(genome_name):
    """Recreate bed files with all important information
    imported file format:
    585	NM_001005484.2	chr1	+	65418	71585	65564	70008	3	65418,65519,69036,	65433,65573,71585,	0	OR4F5	cmpl	cmpl	-1,0,0,

    """
    genome_dir = os.path.join(
            os.path.dirname(__file__), 'localization_data', 'temp')
    sorted_genome_file = os.path.join(genome_dir, f'{genome_name}_sorted.txt')
    gene_localization_file_name = os.path.join(
        os.path.dirname(__file__), 'localization_data', 'genes_coordinates',
        f'{genome_name}.bed')
    chromosome_data = defaultdict(list)
    with open(sorted_genome_file) as f_in:
        with open(gene_localization_file_name, 'w') as f_genes:
            new_data = {}
            for line in f_in:
                elements = line.strip().split()
                chromosome = elements[2]
                if len(chromosome) > 7:
                    continue
                #chr1	11873	14409	NR_046018.2	0	+	14409	14409
                # new:
                # chromosome, position_start, position_end, gene, 0, strand,
                pos_start = int(elements[4])
                pos_end = int(elements[5])
                if elements[12] in new_data:
                    pos_start = min(pos_start, new_data[elements[12]][1])
                    pos_end = max(pos_end, new_data[elements[12]][2])

                new_data[elements[12]] = (
                    chromosome, pos_start, pos_end, elements[12], 0,
                    elements[3])
                # ADD CHROMOSOME DATA (FOR ELEMENTS, UTRS, CODING ETC.)
                # here add chr3 65641 66175 NR_110824.1_utr5_1_0_chr3_65642_f 5utr

                transcription_start = int(elements[4])
                transcription_end = int(elements[5])
                coding_start = int(elements[6])
                coding_end = int(elements[7])

                coding_start_sites = elements[9].split(',')
                coding_end_sites = elements[10].split(',')
                if coding_start_sites[-1] == '':
                    coding_start_sites.pop()
                if coding_end_sites[-1] == '':
                    coding_end_sites.pop()
                coding_start_sites = [int(x) for x in coding_start_sites]
                coding_end_sites = [int(x) for x in coding_end_sites]
                for i in range(len(coding_start_sites)):
                    current_utr_5_start = None
                    current_utr_5_end = None
                    current_utr_3_start = None
                    current_utr_3_end = None
                    current_coding_start = None
                    current_coding_end = None
                    if coding_start <= coding_start_sites[i]:
                        current_coding_start = coding_start_sites[i]

                        if coding_end >= coding_end_sites[i]:
                            current_coding_end = coding_end_sites[i]
                        else:
                            if coding_end > coding_start_sites[i]:
                                current_coding_end = coding_end
                                current_utr_3_start = coding_end + 1
                            else:
                                current_utr_3_start = coding_start_sites[i]
                            current_utr_3_end = coding_end_sites[i]
                    else:
                        current_utr_5_start = coding_start_sites[i]
                        if coding_start < coding_end_sites[i]:
                            current_utr_5_end = coding_start -1
                            current_coding_start = coding_start
                            if coding_end >= coding_end_sites[i]:
                                current_coding_end = coding_end_sites[i]
                            else:
                                current_coding_end = coding_end
                                current_utr_3_start = coding_end + 1
                                current_utr_3_end = coding_end_sites[i]
                        else:
                            current_utr_5_end = coding_end_sites[i]

                    # FIX for UCSC notation - adding +1 to start
                    if current_utr_5_start and current_utr_5_end:
                        chromosome_data[chromosome].append((
                            current_utr_5_start+1, current_utr_5_end,
                            elements[12], 'utr')) #elements[1] if NM_015102
                    if current_coding_start and current_coding_end:
                        chromosome_data[chromosome].append((
                            current_coding_start+1, current_coding_end,
                            elements[12], 'coding_exons'))#elements[1] if NM_015102
                    if current_utr_3_start and current_utr_3_end:
                        chromosome_data[chromosome].append((
                            current_utr_3_start+1, current_utr_3_end,
                            elements[12], 'utr'))#elements[1] if NM_015102
                    # introns
                    if i + 1 != len(coding_end_sites):
                        chromosome_data[chromosome].append((
                            coding_end_sites[i]+1, coding_start_sites[i + 1],
                            elements[12], 'introns'))#elements[1] if NM_015102

            sorted_new_data = sorted(
                list(new_data.values()),
                key=lambda x: (not x[0][3:].isdigit(),
                               len(x[0]), x[0], x[1] if x[-1] == '+' else x[2],
                               x[2] if x[-1] == '+' else x[1]))
            for entry in sorted_new_data:
                f_genes.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(*entry))

    # write data
    for chromosome, values in chromosome_data.items():
        filename = os.path.join(
            os.path.dirname(__file__), 'localization_data', 'elements_coordinates',
            f"{genome_name}_{chromosome}.bed")
        values = list(set(values))
        sorted_values = sorted(values, key=lambda x: (x[0], x[1], x[2]))
        with open(filename, 'w') as f:
            for entry in sorted_values:
                f.write(f'{chromosome} {entry[0]} {entry[1]} '
                        f'{entry[2]} {entry[3]}\n')
    # cleanup
    os.remove(sorted_genome_file)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('genomes', metavar='N', type=str, nargs='+',
                    help='genome names')
    args = parser.parse_args()
    create_directories()
    for genome in args.genomes:
        download_single_genome_data(genome)
        create_gene_and_elements_files(genome)


if __name__ == '__main__':
    main()














