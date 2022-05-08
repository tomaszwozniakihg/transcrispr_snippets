These are selected parts of code used in transCRISPR project. 
prepare_localization.py will create required by localization.py custom .bed files
usage:
python3.8 prepare_localization.py hg19 hg38

localization.py will find localization of given motifs
usage:
python3.8 localization.py -i motifs.csv -o o.txt -g hg19

find_motif_in_sequence.py is code that can be used to find specific motif and guide position
extended sequence is required if you want to find guides for motifs that are at the begining or end of the sequence
usage:
python3.8 find_motif_in_sequence.py -variant 1 -chromosome chr2 -sequence sequence.txt -motifs motifs.txt -extended extended.txt -binding_site site1 -binding_start 264308 -binding_end 264501

utils.py and exceptions.py are files with code parts that one might find usefull

biopython_mocks contain tools to transform motifs from biopython into IUPAC code
