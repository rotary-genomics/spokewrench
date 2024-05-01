#!/usr/bin/env python
# external.py
# Code for handling external CLI applications used by rotary-utils (e.g. flye)
# Copyright Jackson M. Tsuji and Lee H. Bergstrand 2024
import logging
import os
import shlex
import shutil
import subprocess

import pandas as pd

from rotary_utils.utils import set_write_mode

logger = logging.getLogger(__name__)


def check_dependency(dependency_name: str):
    """
    Checks if a required shell dependency is present and tries to get its version.

    :param dependency_name: name of the dependency
    :return: tuple of the path to the dependency and dependency version
    """

    dependency_path = shutil.which(dependency_name)
    if dependency_path is None:
        error = FileNotFoundError(f'Dependency not found: {dependency_name}')
        logger.error(error)
        raise error

    dependency_version = get_dependency_version(dependency_name)

    output_tuple = (dependency_path, dependency_version)
    return output_tuple


def check_dependencies(dependency_names: list):
    """
    For each provided dependency name, checks if the dependency exists and gets the path and version.

    :param dependency_names: a list of names of dependencies to check
    :return: dictionary with dependency names as key and a tuple of dependency paths and versions as values
    """

    path_and_version_tuples = []
    for dependency_name in dependency_names:
        path_and_version_tuple = check_dependency(dependency_name)
        path_and_version_tuples.append(path_and_version_tuple)

        dependency_path, dependency_version = path_and_version_tuple
        logger.debug(f'{dependency_name}: version {dependency_version}: {dependency_path}')

    dependency_dict = dict(zip(dependency_names, path_and_version_tuples))

    return dependency_dict


def get_dependency_version(dependency_name: str, log: bool = False):
    """
    Tries to get the version of a dependency based on the name of the dependency.
    :param dependency_name: name of the dependency
    :param log: is True, print a log of the shell commands used
    :return: version of the dependency. 'unknown' if no version can be parsed
    """
    version_commands = {
        'flye': ['flye', '--version'],
        'minimap2': ['minimap2', '--version'],
        'samtools': ['samtools', 'version'],
        'circlator': ['circlator', 'version'],
    }

    if dependency_name in version_commands:
        stdout = run_pipeline_subcommand(version_commands[dependency_name], stdout=subprocess.PIPE,
                                         text=True, log=log).stdout

        # Special parsing is needed for samtools output.
        if dependency_name == 'samtools':
            dependency_version = stdout.splitlines()[0].split(' ')[1]
        else:
            dependency_version = stdout.splitlines()[0]
    else:
        dependency_version = 'unknown'
    return dependency_version


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


def map_long_reads(contig_filepath: str, long_read_filepath: str, output_bam_filepath: str, log_filepath: str,
                   append_log: bool = True, threads: int = 1, threads_mem_mb: float = 1):
    """
    Maps long reads (via minimap2) to contigs and sorts/indexes the resulting BAM file. output_bam_filepath and
    log_filepath are saved to disk.

    :param contig_filepath: path to the FastA file containing the reference contigs
    :param long_read_filepath: path to the FastQ file containing long reads to map (compressed is OK)
    :param output_bam_filepath: path to the BAM file to be saved
    :param log_filepath: path to the log file to be saved
    :param append_log: whether the log should append to an existing file (True) or overwrite an existing file (False);
                       this setting is only relevant if the log file at log_filepath already exists
    :param threads: number of threads to use for read mapping
    :param threads_mem_mb: memory in MB per thread (to use for samtools); must be an integer
    """

    write_mode = set_write_mode(append_log)

    with open(log_filepath, write_mode) as logfile_handle:
        with open(output_bam_filepath, 'w') as bam_handle:
            # TODO - add support for different flags like -ax for pacbio
            minimap_args = ['minimap2', '-t', str(threads), '-ax', 'map-ont', contig_filepath, long_read_filepath]
            minimap = run_pipeline_subcommand(command_args=minimap_args, stdout=subprocess.PIPE, stderr=logfile_handle)

            samtools_view_args = ['samtools', 'view', '-b', '-@', str(threads)]
            samtools_view = run_pipeline_subcommand(command_args=samtools_view_args, stdin=minimap,
                                                    stdout=subprocess.PIPE, stderr=logfile_handle)

            samtools_sort_args = ['samtools', 'sort', '-@', str(threads), '-m', f'{threads_mem_mb}M']
            run_pipeline_subcommand(command_args=samtools_sort_args, stdin=samtools_view, stdout=bam_handle,
                                    stderr=logfile_handle)

        samtools_index_args = ['samtools', 'index', '-@', str(threads), output_bam_filepath]
        run_pipeline_subcommand(command_args=samtools_index_args, stderr=logfile_handle)

    logger.debug('Read mapping finished')


def subset_reads_from_bam(bam_filepath: str, bed_filepath: str, subset_fastq_filepath: str, log_filepath: str,
                          append_log: bool = True, threads: int = 1):
    """
    Subsets reads from a BAM file that were mapped to regions defined in a BED file; saves reads to a FastQ file, at
    filepath subset_fastq_filepath.

    :param bam_filepath: path to a BAM file containing reads mapped to a reference; BAM needs to be sorted and indexed
    :param bed_filepath: path to a BED file containing the regions of reference contigs to subset reads for
    :param subset_fastq_filepath: path to the FastQ file to be saved (.fastq.gz extension saves as Gzipped FastQ)
    :param log_filepath: path to the log file to be saved
    :param append_log: whether the log should append to an existing file (True) or overwrite an existing file (False);
                       this setting is only relevant if the log file at log_filepath already exists
    :param threads: number of threads to use for read mapping
    """

    write_mode = set_write_mode(append_log)

    with open(log_filepath, write_mode) as logfile_handle:
        # TODO - consider adding option to split long reads in half if they go around a short circular contig,
        #  like in circlator
        samtools_view_args = ['samtools', 'view', '-@', str(threads), '-L', bed_filepath, '-b', bam_filepath]
        samtools_view = run_pipeline_subcommand(command_args=samtools_view_args, stdout=subprocess.PIPE,
                                                stderr=logfile_handle)

        samtools_fastq_args = ['samtools', 'fastq', '-0', subset_fastq_filepath, '-n', '-@', str(threads)]
        run_pipeline_subcommand(command_args=samtools_fastq_args, stdin=samtools_view, stderr=logfile_handle)


def run_flye(fastq_filepath: str, flye_outdir: str, flye_read_mode: str, flye_read_error: float, log_filepath: str,
             append_log: bool = True, threads: int = 1):
    """
    Runs Flye to assemble the reads in the input FastQ file. This function allows Flye to fail without raising an error.

    :param fastq_filepath: path to a FastQ file containing the input reads (gzipped is OK)
    :param flye_outdir: directory to save Flye output to
    :param flye_read_mode: type of long reads, either 'nano-raw' or 'nano-hq'
    :param flye_read_error: expected error rate of reads as a proportion; specify 0 to use default Flye settings
    :param log_filepath: path to the log file to be saved
    :param append_log: whether the log should append to an existing file (True) or overwrite an existing file (False);
                       this setting is only relevant if the log file at log_filepath already exists
    :param threads: number of threads to use for read mapping
    :return: return code of Flye (0 if it finished successfully).
    """

    # TODO - add support for PacBio reads
    if (flye_read_mode != 'nano-raw') & (flye_read_mode != 'nano-hq'):
        error = ValueError(f'flye_read_mode must be "nano-raw" or "nano-hq"; you provided {flye_read_mode}')
        logger.error(error)
        raise error

    flye_args = ['flye', f'--{flye_read_mode}', fastq_filepath, '-o', flye_outdir, '-t', str(threads)]

    if flye_read_error != 0:
        flye_args.append('--read_error')
        flye_args.append(flye_read_error)

    write_mode = set_write_mode(append_log)
    with open(log_filepath, write_mode) as logfile_handle:
        flye_run = run_pipeline_subcommand(command_args=flye_args, check=False, stderr=logfile_handle)

    if flye_run.returncode != 0:
        logger.warning(f'Flye did not finish successfully; see log for details at "{log_filepath}"')

    return flye_run.returncode


def run_circlator_merge(circular_contig_filepath: str, patch_contig_filepath: str, merge_outdir: str,
                        circlator_min_id: float, circlator_min_length: int, circlator_ref_end: int,
                        circlator_reassemble_end: int, log_filepath: str, append_log: bool = True):
    """
    Runs the 'circlator merge' module to stitch a gap-spanning contig onto the ends of a circular contig to confirm and
    repair the circularization of the contig. 'circlator merge' output is saved to disk at merge_outdir

    :param circular_contig_filepath: path to a FastA file containing the original (non-stitched) circular contig
    :param patch_contig_filepath: path to a FastA file containing a linear contig that should span the 'ends' of the
                                  circular contig
    :param merge_outdir: directory to save circlator merge output to
    :param circlator_min_id: Percent identity threshold for circlator merge to stitch the contigs
    :param circlator_min_length: Minimum required overlap (bp) between the circular contig and the patch contig
    :param circlator_ref_end: Minimum distance (bp) between end of circular contig and the nucmer hit
    :param circlator_reassemble_end: Minimum distance (bp) between end of patch contig and the nucmer hit
    :param log_filepath: path to the log file to be saved
    :param append_log: whether the log should append to an existing file (True) or overwrite an existing file (False);
                       this setting is only relevant if the log file at log_filepath already exists
    """

    os.makedirs(merge_outdir, exist_ok=True)

    write_mode = set_write_mode(append_log)
    with open(log_filepath, write_mode) as logfile_handle:
        circlator_merge_args = ['circlator', 'merge', '--verbose', '--min_id', str(circlator_min_id),
                                '--min_length', str(circlator_min_length), '--ref_end', str(circlator_ref_end),
                                '--reassemble_end', str(circlator_reassemble_end), circular_contig_filepath,
                                patch_contig_filepath, os.path.join(merge_outdir, 'merge')]
        run_pipeline_subcommand(command_args=circlator_merge_args, stdout=logfile_handle, stderr=subprocess.STDOUT)


def check_circlator_success(circlator_logfile: str):
    """
    Checks the circlator log file to see if contig stitching was successful.

    :param circlator_logfile: path to the circlator log file
    :return: Boolean of whether the contigs were successfully stitched or not
    """

    logger.debug('Loading circlator logfile')
    circlator_info = pd.read_csv(circlator_logfile, sep='\t')[['#Contig', 'circularised']]

    # TODO: I think I might be able to safely just look for if the sum is 1 (>1 might mean something odd happened)
    # If a row is '1', it means it was stitched properly, but if '0', it means it was not stitched.
    # So if all rows are 1 (i.e., the sum of rows / # of rows is 1), it means everything was stitched properly.
    if circlator_info['circularised'].sum() / circlator_info.shape[0] == 1:
        logger.debug('Everything is stitched.')
        result = True
    elif circlator_info['circularised'].sum() >= 0:
        logger.debug('Not everything is stitched.')
        result = False
    else:
        error = RuntimeError(f'Could not understand the circularization status output by circlator: '
                             f'{circlator_info.shape[0]} rows, but the sum of circularization status of these rows is '
                             f'{circlator_info["circularised"].sum()} (should be >=0, because each row should have a'
                             f'status of either 0 or 1).')
        logger.error(error)
        raise error

    return result
