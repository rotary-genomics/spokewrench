#!/usr/bin/env python
# utils.py
# Utility function within rotary-utils
# Copyright Jackson M. Tsuji and Lee H. Bergstrand 2024

import logging
import shlex
import shutil
import subprocess

logger = logging.getLogger(__name__)


def check_dependency(dependency_name: str):
    """
    Checks if a required shell dependency is present
    :param dependency_name: name of the dependency
    :return: path to the dependency
    """

    dependency_path = shutil.which(dependency_name)

    if dependency_path is None:

        logger.error(f'Dependency not found: {dependency_name}')
        raise RuntimeError

    return dependency_path


def set_write_mode(append_log: bool):
    """
    Converts the boolean append_log to 'w' or 'a' write modes
    :param append_log: boolean of whether to append to an existing log file (True) or to overwrite an existing log
                       file (False)
    :return: string of either 'a' (append mode) or 'w' (write mode)
    """

    if append_log is True:
        write_mode = 'a'

    elif append_log is False:
        write_mode = 'w'

    else:

        logger.error(f'append_log should be a boolean True or False; you provided {append_log}')
        raise ValueError

    return write_mode


def run_pipeline_subcommand(command_args, stdin=None, stdout=None, stderr=None, check=True):
    """
    Wrapper function for running subcommands.

    :param command_args: The command line arguments of the subcommand (e.g., ['samtools', '-h'])
    :param stdin: A subprocess.PIPE or None if stdin is not to be used.
    :param stdout: Where to send stdout or None if stdout is not to be used.
    :param stderr: Where to send stderr or None if stderr is not to be used.
    :param check: Cause a runtime error if the subcommand fails.
    :return: The output of the subcommand.
    """
    logger.debug(shlex.join(command_args))

    if stdin:
        stdin = stdin.stdout

    return subprocess.run(command_args, check=check, input=stdin, stdout=stdout, stderr=stderr)


def load_fasta_sequences(fasta_filepath: str, maximum_allowed_sequences: int = 0):
    """
    Loads FastA sequence records as BioPython SeqRecord objects.
    :param fasta_filepath: Path to the FastA file (unzipped) to load
    :param maximum_allowed_sequences: Maximum number of sequences allowed in the input file to load without error.
                                      If 0, then no sequence number limit is applied.
    :return: list of SeqRecord objects (one object per sequence loaded)
    """

    record_count = 0
    sequence_records = []

    with open(fasta_filepath) as fasta_handle:
        for record in SeqIO.parse(fasta_handle, 'fasta'):

            if (record_count < maximum_allowed_sequences) | (maximum_allowed_sequences == 0):
                sequence_records.append(record)

            elif record_count >= maximum_allowed_sequences:
                logger.error(f'Input file {fasta_filepath} contains more than the maximum record limit of'
                             f'{maximum_allowed_sequences} sequences.')
                raise RuntimeError

            record_count = record_count + 1

    return sequence_records
