class ImproperNucleotideException(Exception):
    """One of the nucleotides in the sequence is improper"""


class ImproperCoordinateFileException(Exception):
    """File with coordinates is formatted incorrectly use:
    'chrX start end' format"""


class NoSequenceException(Exception):
    """Exception raised when sequence cannot be found"""


