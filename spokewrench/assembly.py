#!/usr/bin/env python
# assembly.py
# Module containing functions related to processing Flye-based sequence assembly
# Copyright Jackson M. Tsuji and Lee H. Bergstrand 2024

import logging
import pandas as pd

logger = logging.getLogger(__name__)


class AssemblyInfo:
    """
    A class representing file paths containing info about the input assembly.
    """

    def __init__(self, assembly_fasta_filepath: str, assembly_info_filepath: str, assembly_info_type: str):
        """
        Instantiate an AssemblyInfo object.

        :param assembly_fasta_filepath: path to the assembled contigs output by Flye, FastA format
        :param assembly_info_filepath: path to assembly_info.txt output by Flye
        :param assembly_info_type: 'flye' type for assembly_info.txt from flye; 'custom' for custom format (see
                                   parse_cli for custom format details)
        """
        self.assembly_fasta_filepath = assembly_fasta_filepath
        self.assembly_info_filepath = assembly_info_filepath
        self.assembly_info_type = assembly_info_type


def load_custom_assembly_info_file(assembly_info_filepath: str):
    """
    Load and check a custom-format assembly info file.

    :param assembly_info_filepath: path to the custom assembly info file.
    :return: pandas DataFrame of the assembly info file
    """

    assembly_info = pd.read_csv(assembly_info_filepath, sep='\t', header=None)
    assembly_info.columns = ['#seq_name', 'status']

    expected_statuses = {'circular', 'linear'}
    actual_statuses = set(assembly_info['status'])
    unexpected_statuses = actual_statuses.difference(expected_statuses)

    if len(unexpected_statuses) != 0:
        logger.warning(f'Some entries in the assembly info file had unexpected contig status names, i.e.: '
                       f'{", ".join(unexpected_statuses)}')
        logger.warning('These entries will be treated as linear contigs... they will not be rotated and will be '
                       'returned as-is at the end of the script. Please make sure you did not make a typo or '
                       'include a header for your custom assembly info file.')

    return assembly_info


def parse_assembly_info_file(assembly_info_filepath: str, info_type: str):
    """
    Extract circular and linear contig names from a Flye (or custom format) assembly info file.

    :param assembly_info_filepath: path to assembly_info.txt output by Flye or a custom assembly info file.
    :param info_type: whether the info file is in 'flye' format or is a 'custom' format.
                      'flye' format refers to the 'assembly_info.txt' format output by Flye after a successful assembly.
                      'custom' info files are tab-separated, have no headers, and have two columns: contig name and
                      contig type, either 'circular' or 'linear'.
    :return: tuple containing a list of circular contig names (entry 0) and a list of linear contig names (entry 1).
    """

    if info_type == 'flye':
        logger.debug('Loading Flye-format assembly info file')
        assembly_info = pd.read_csv(assembly_info_filepath, sep='\t')
        circular_contigs = assembly_info[assembly_info['circ.'] == 'Y']
        linear_contigs = assembly_info[assembly_info['circ.'] == 'N']
    elif info_type == 'custom':
        logger.debug('Loading custom format assembly info file')
        assembly_info = load_custom_assembly_info_file(assembly_info_filepath)
        circular_contigs = assembly_info[assembly_info['status'] == 'circular']
        linear_contigs = assembly_info[assembly_info['status'] != 'circular']
    else:
        error = ValueError(f'Assembly info file format should be "flye" or "custom"; you provided {info_type}')
        logger.error(error)
        raise error

    # Check for duplicate sequence IDs
    if assembly_info['#seq_name'].drop_duplicates().shape[0] < assembly_info['#seq_name'].shape[0]:
        duplicated_ids = set(assembly_info['#seq_name'][assembly_info['#seq_name'].duplicated() == True])
        error = RuntimeError(f'Some sequence IDs are duplicated in the input assembly info file: '
                             f'{", ".join(duplicated_ids)}')
        logger.error(error)
        raise error

    contig_summary = (list(circular_contigs['#seq_name']), list(linear_contigs['#seq_name']))
    logger.debug(f'Found {len(contig_summary[0])} circular sequences and {len(contig_summary[1])} linear sequences')

    return contig_summary
