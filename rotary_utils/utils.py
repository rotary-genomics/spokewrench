#!/usr/bin/env python
# utils.py
# Utility function within rotary-utils
# Copyright Jackson M. Tsuji and Lee H. Bergstrand 2024
import logging
import os
import sys

import pandas as pd
from Bio import SeqIO

logger = logging.getLogger(__name__)


def check_output_file(output_filepath: str, overwrite: bool = False):
    """
    Checks if OK to create an output file. Raises an error if the output file already exists (unless overwrite=True).

    :param output_filepath: path to the desired output file
    :param overwrite: if True, the keep going with a warning if the output file already exists
    """

    output_file_exists = os.path.isfile(output_filepath)

    if output_file_exists is True:
        if overwrite is False:
            logger.error(f'Output file already exists: "{output_filepath}". Will not continue. Set the '
                         f'--overwrite flag at your own risk if you want to overwrite existing files.')
            sys.exit(1)
        elif overwrite is True:
            logger.warning(f'Output file already exists: "{output_filepath}". File will be overwritten.')
        else:
            raise ValueError(f'overwrite must be True or False, but you provided "{overwrite}"')


def set_up_output_directory(output_directory_filepath: str, overwrite: bool = False):
    """
    Creates an output directory. Raises an error if a directory already exists (unless overwrite=True).

    :param output_directory_filepath: path to the desired output directory
    :param overwrite: if True, then keep going with a warning if the output directory already exists
    """

    output_dir_exists = os.path.isdir(output_directory_filepath)

    if output_dir_exists is True:
        if overwrite is False:
            logger.error(f'Output directory already exists: "{output_directory_filepath}". Will not continue. Set the '
                         f'--overwrite flag at your own risk if you want to use an existing directory.')
            sys.exit(1)
        elif overwrite is True:
            logger.warning(f'Output directory already exists: "{output_directory_filepath}". Files may be overwritten.')
        else:
            raise ValueError(f'overwrite must be True or False, but you provided "{overwrite}"')

    os.makedirs(output_directory_filepath, exist_ok=True)


def set_write_mode(append_log: bool):
    """
    Converts the boolean append_log to 'w' or 'a' write modes.

    :param append_log: boolean of whether to append to an existing log file (True) or to overwrite an existing log
                       file (False)
    :return: string of either 'a' (append mode) or 'w' (write mode)
    """

    if append_log is True:
        write_mode = 'a'
    elif append_log is False:
        write_mode = 'w'
    else:
        error = ValueError(f'append_log should be a boolean True or False; you provided {append_log}')
        logger.error(error)
        raise error

    return write_mode


def load_fasta_sequences(fasta_filepath: str):
    """
    Loads an input FastA file as a generator.

    :param fasta_filepath: Path to the FastA file (unzipped) to load
    :return: generator of a SeqRecord object for the loaded sequences
    """

    with open(fasta_filepath) as fasta_handle:
        for record in SeqIO.parse(fasta_handle, 'fasta'):
            yield record


def subset_sequences(input_fasta_filepath: str, subset_sequence_ids: list):
    """
    Given an input FastA file, subsets the file to the provided sequence IDs.

    :param input_fasta_filepath: Path to the input FastA file
    :param subset_sequence_ids: list of the names of sequences to keep. If any names in the list are not in the
                                input file, the function will not return anything for those names. The function will
                                raise an error if any duplicate sequence names are detected in the input FastA file.
    :return: generator of a SeqRecord object for the subset sequences
    """

    try:
        sequence_names = []
        for record in load_fasta_sequences(input_fasta_filepath):
            sequence_names.append(record.name)

            if record.name in subset_sequence_ids:
                yield record
    finally:
        # Raise an error if there are duplicate sequence names
        if len(set(sequence_names)) < len(sequence_names):
            sequence_names_series = pd.Series(sequence_names)
            duplicates_names = set(sequence_names_series[sequence_names_series.duplicated() == True])

            error = RuntimeError(f'Duplicate sequence IDs were detected in the input FastA file '
                                 f'"{input_fasta_filepath}": {", ".join(duplicates_names)}')
            logger.error(error)
            raise error
