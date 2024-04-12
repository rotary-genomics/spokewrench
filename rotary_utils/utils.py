#!/usr/bin/env python
# utils.py
# Utility function within rotary-utils
# Copyright Jackson M. Tsuji and Lee H. Bergstrand 2024

import os
import sys
import logging
import shlex
import shutil
import subprocess

from Bio import SeqIO

logger = logging.getLogger(__name__)


def check_output_file(output_filepath: str, overwrite: bool = False):
    """
    Checks if OK to create an output file. Raises an error if the output file already exists (unless overwrite=True)
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
            raise ValueError


def set_up_output_directory(output_directory_filepath: str, overwrite: bool = False):
    """
    Creates an output directory. Raises an error if a directory already exists (unless overwrite=True)
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
            raise ValueError

    os.makedirs(output_directory_filepath, exist_ok=True)


def get_dependency_version(dependency_name: str, log: bool = False):
    """
    Tries to get the version of a dependency based on the name of the dependency
    :param dependency_name: name of the dependency
    :param log: is True, print a log of the shell commands used
    :return: version of the dependency. 'unknown' if no version can be parsed
    """

    if dependency_name == 'flye':
        stdout = run_pipeline_subcommand(['flye', '--version'], stdout=subprocess.PIPE, text=True,
                                         log=log).stdout
        dependency_version = stdout.splitlines()[0]

    elif dependency_name == 'minimap2':
        stdout = run_pipeline_subcommand(['minimap2', '--version'], stdout=subprocess.PIPE, text=True,
                                         log=log).stdout
        dependency_version = stdout.splitlines()[0]

    elif dependency_name == 'samtools':
        stdout = run_pipeline_subcommand(['samtools', 'version'], stdout=subprocess.PIPE, text=True,
                                         log=log).stdout
        dependency_version = stdout.splitlines()[0].split(' ')[1]

    elif dependency_name == 'circlator':
        stdout = run_pipeline_subcommand(['circlator', 'version'], stdout=subprocess.PIPE, text=True,
                                         log=log).stdout
        dependency_version = stdout.splitlines()[0]

    else:
        dependency_version = 'unknown'

    return dependency_version


def check_dependency(dependency_name: str):
    """
    Checks if a required shell dependency is present and tries to get its version
    :param dependency_name: name of the dependency
    :return: tuple: path to the dependency, dependency version (if available)
    """

    dependency_path = shutil.which(dependency_name)

    if dependency_path is None:

        logger.error(f'Dependency not found: {dependency_name}')
        raise RuntimeError

    dependency_version = get_dependency_version(dependency_name)

    output_tuple = (dependency_path, dependency_version)
    return output_tuple


def check_dependencies(dependency_names: list):
    """
    For each provided dependency name, checks if the dependency exists and gets the path and version.
    :param dependency_names: a list of names of dependencies to check
    :return: a dictionary with dependency names as key and a tuple of dependency paths and versions as values
    """

    path_and_version_tuples = []
    for dependency_name in dependency_names:
        path_and_version_tuple = check_dependency(dependency_name)
        path_and_version_tuples.append(path_and_version_tuple)

        dependency_path, dependency_version = path_and_version_tuple
        logger.debug(f'{dependency_name}: version {dependency_version}: {dependency_path}')

    dependency_dict = dict(zip(dependency_names, path_and_version_tuples))

    return dependency_dict


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


def run_pipeline_subcommand(command_args, stdin=None, stdout=None, stderr=None, check=True, text=None, log=True):
    """
    Wrapper function for running subcommands.

    :param command_args: The command line arguments of the subcommand (e.g., ['samtools', '-h'])
    :param stdin: A subprocess.PIPE or None if stdin is not to be used.
    :param stdout: Where to send stdout or None if stdout is not to be used.
    :param stderr: Where to send stderr or None if stderr is not to be used.
    :param check: Cause a runtime error if the subcommand fails.
    :param text: If True, open outputs in text mode (rather than binary)
    :param log: If True, write the shell command to logger in debug mode
    :return: The output of the subcommand.
    """
    if log is True:
        logger.debug(shlex.join(command_args))

    if stdin:
        stdin = stdin.stdout

    return subprocess.run(command_args, check=check, input=stdin, stdout=stdout, stderr=stderr, text=text)


def load_fasta_sequences(fasta_filepath: str):
    """
    Loads an input FastA file as a generator

    :param fasta_filepath: Path to the FastA file (unzipped) to load
    :return: generator of a SeqRecord object for the loaded sequences
    """

    with open(fasta_filepath) as fasta_handle:
        for record in SeqIO.parse(fasta_handle, 'fasta'):
            yield record


def load_all_fasta_sequences(fasta_filepath: str, maximum_allowed_sequences: int = 0):
    """
    Loads FastA sequence records as BioPython SeqRecord objects. Rather than yielding a generator as in
    load_fasta_sequences(), this function returns a list of all loaded sequence records in RAM.

    :param fasta_filepath: Path to the FastA file (unzipped) to load
    :param maximum_allowed_sequences: Maximum number of sequences allowed in the input file to load without error.
                                      If 0, then no sequence number limit is applied.
    :return: list of SeqRecord objects (one object per sequence loaded)
    """

    record_count = 0
    sequence_records = []

    for record in load_fasta_sequences(fasta_filepath):

        if (record_count < maximum_allowed_sequences) | (maximum_allowed_sequences == 0):
            sequence_records.append(record)

        elif record_count >= maximum_allowed_sequences:
            logger.error(f'Input file {fasta_filepath} contains more than the maximum record limit of'
                         f'{maximum_allowed_sequences} sequences.')
            raise RuntimeError

        record_count = record_count + 1

    return sequence_records
