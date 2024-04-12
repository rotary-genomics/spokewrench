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

logger = logging.getLogger(__name__)


def check_log_file(log_filepath: str, overwrite: bool = False):
    """
    Checks if OK to create a log file. Raises an error if the log file already exists (unless overwrite=True)
    :param log_filepath: path to the desired log file
    :param overwrite: if True, the keep going with a warning if the log file already exists
    """

    log_file_exists = os.path.isfile(log_filepath)

    if log_file_exists is True:
        if overwrite is False:
            logger.error(f'Log file already exists: "{log_filepath}". Will not continue. Set the '
                         f'--overwrite flag at your own risk if you want to overwrite existing files.')
            sys.exit(1)

        elif overwrite is True:
            logger.warning(f'Log file already exists: "{log_filepath}". File will be overwritten.')

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


def check_dependencies(dependency_names: list):
    """
    For each provided dependency name, checks if the dependency exists and gets the path.
    :param dependency_names: a list of names of dependencies to check
    :return: a dictionary of dependency names and dependency paths
    """

    dependency_paths = []
    for dependency_name in dependency_names:
        dependency_path = check_dependency(dependency_name)
        dependency_paths.append(dependency_path)
        logger.debug(f'{dependency_name}: {dependency_path}')

    dependency_dict = dict(zip(dependency_names, dependency_paths))

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
