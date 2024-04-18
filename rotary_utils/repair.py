#!/usr/bin/env python
# repair.py
# Fixes ends of circular contigs produced by Flye
# Copyright Jackson M. Tsuji and Lee H. Bergstrand 2024

import argparse
import logging
import os
import shutil
import sys

import pandas as pd
from Bio import SeqIO

from rotary_utils.external import map_long_reads, subset_reads_from_bam, run_flye, run_circlator_merge, \
    check_circlator_success, check_dependencies
from rotary_utils.rotate import rotate_sequences_wf
from rotary_utils.utils import load_fasta_sequences

# GLOBAL VARIABLES
DEPENDENCY_NAMES = ['flye', 'minimap2', 'samtools', 'circlator']

logger = logging.getLogger(__name__)


def main(args):
    """
    Runs the end repair workflow.

    :param args: arguments parsed by parse_cli()
    """

    # Parse some command line inputs further
    assembly_info_type = 'custom' if args.custom_assembly_info_file is True else 'flye'
    # TODO - improve error handling if a string is provided instead of a real length
    length_thresholds = [int(x) for x in args.length_thresholds.split(',')]
    cli_tool_settings_dict = {'flye_read_mode': args.flye_read_mode,
                              'flye_read_error': args.flye_read_error,
                              'circlator_min_id': args.circlator_min_id,
                              'circlator_min_length': args.circlator_min_length,
                              'circlator_ref_end': args.circlator_ref_end,
                              'circlator_reassemble_end': args.circlator_reassemble_end}
    threads_mem_mb = int(args.threads_mem * 1024)

    # Startup messages
    logger.info('Running ' + os.path.basename(sys.argv[0]))
    logger.debug('### DEPENDENCIES ###')
    check_dependencies(dependency_names=DEPENDENCY_NAMES)
    logger.debug('### SETTINGS ###')
    logger.debug(f'Long read filepath: {args.long_read_filepath}')
    logger.debug(f'Assembly FastA filepath: {args.assembly_fasta_filepath}')
    logger.debug(f'Assembly info filepath: {args.assembly_info_filepath}')
    logger.debug(f'Use a custom assembly info file: {args.custom_assembly_info_file}')
    logger.debug(f'Output directory: {args.output_dir}')
    logger.debug(f'Overwrite output dir if needed?: {args.overwrite}')
    logger.debug(f'Length thresholds to test (bp): {length_thresholds}')
    logger.debug(f'Keep going if some contigs cannot be re-circularized?: {args.keep_going_with_failed_contigs}')
    logger.debug(f'Flye read mode: {cli_tool_settings_dict["flye_read_mode"]}')
    logger.debug(f'Flye read error: {cli_tool_settings_dict["flye_read_error"]}')
    logger.debug(f'Circlator min. ID: {cli_tool_settings_dict["circlator_min_id"]}')
    logger.debug(f'Circlator min. length: {cli_tool_settings_dict["circlator_min_length"]}')
    logger.debug(f'Circlator ref. end: {cli_tool_settings_dict["circlator_ref_end"]}')
    logger.debug(f'Circlator reassembly end: {cli_tool_settings_dict["circlator_reassemble_end"]}')
    logger.debug(f'Custom log file path: {args.logfile}')
    logger.debug(f'Threads: {args.threads}')
    logger.debug(f'Memory per thread (GB = MB): {args.threads_mem} = {threads_mem_mb}')
    logger.debug(f'Verbose logging: {args.verbose}')
    logger.debug('################')

    assembly_info = AssemblyInfo(args.assembly_fasta_filepath, args.assembly_info_filepath, assembly_info_type)

    run_end_repair(args.long_read_filepath, assembly_info, args.output_dir, length_thresholds,
                   args.keep_going_with_failed_contigs, cli_tool_settings_dict, args.threads, threads_mem_mb)

    logger.info(os.path.basename(sys.argv[0]) + ': done.')


class AssemblyInfo:
    """
    A class representing files paths containing info about the input assembly.
    """

    def __init__(self, assembly_fasta_filepath, assembly_info_filepath, assembly_info_type):
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


class RepairPaths:
    """
    A class containing file paths that are used during the assembly repair process.
    """

    def __init__(self, output_dir):
        """
        Instantiate a RepairPaths object.

        :param output_dir: output directory for the assembly repair process.
        """
        """
        Description of attributes:
            linking_outdir_base: Temporary directory where re-assembly and stitching files will be saved.
            circular_contig_tmp_fasta: Temporary FastA file containing only circular contigs to be stitched.
            bam_filepath: path to a BAM file with mapping information of long reads to the contigs (it is OK if this
                          file also contains mappings to other contigs outside those in the input file).
            circlator_logs: path where a summary of log files from circlator will be saved.
            end_repaired_contigs_filepath: path where the end-repaired contigs will be saved.
            end_repair_report_filepath: path where the end repair report will be saved.
            verbose_logfile: path to a logfile where shell script logs will be saved.
        """
        self.linking_outdir_base = os.path.join(output_dir, 'contigs')
        self.circular_contig_tmp_fasta = os.path.join(self.linking_outdir_base, 'circular_input.fasta')
        self.bam_filepath = os.path.join(output_dir, 'long_read.bam')
        self.circlator_logs = os.path.join(output_dir, 'circlator_logs')
        self.end_repaired_contigs_filepath = os.path.join(output_dir, 'repaired.fasta')
        self.end_repair_report_filepath = os.path.join(output_dir, 'repaired_info.tsv')
        self.verbose_logfile = os.path.join(output_dir, 'verbose.log')


class ContigInfo:
    """
    A class representing file paths containing info about the contigs of an assembly.
    """
    def __init__(self, circular_contig_names, failed_contig_names, linear_contig_names):
        """
        Instantiate a ContigInfo object.

        :param circular_contig_names: list of names of circular contigs in the assembly info file that could be
                                      properly repaired.
        :param failed_contig_names: list of names of circular contigs in the assembly info file that could not be
                                    properly repaired.
        :param linear_contig_names: list of names of non-circular contigs in the assembly info file.
        """
        self.circular_contig_names = circular_contig_names
        self.repaired_contig_names = list(set(circular_contig_names).difference(set(failed_contig_names)))
        self.failed_contig_names = failed_contig_names
        self.linear_contig_names = linear_contig_names


def parse_assembly_info_file(assembly_info_filepath: str, info_type: str):
    """
    List circular and linear contigs from a Flye (or custom format) assembly info file.

    :param assembly_info_filepath: path to assembly_info.txt output by Flye
    :param info_type: whether the info file is in 'flye' format or is a 'custom' format.
                      'flye' format refers to the 'assembly_info.txt' format output by Flye after a successful assembly.
                      'custom' info files are tab-separated, have no headers, and have two columns: contig name and
                      contig type, either 'circular' or 'linear'.
    :return: tuple containing a list of circular contig names (entry 0) and a list of linear contig names (entry 1)
    """

    if info_type == 'flye':
        logger.debug('Loading Flye-format assembly info file')
        assembly_info = pd.read_csv(assembly_info_filepath, sep='\t')
        circular_contigs = assembly_info[assembly_info['circ.'] == 'Y'][['#seq_name', 'length', 'circ.']]
        linear_contigs = assembly_info[assembly_info['circ.'] == 'N'][['#seq_name', 'length', 'circ.']]
    elif info_type == 'custom':
        logger.debug('Loading custom format assembly info file')
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


def subset_sequences(input_fasta_filepath: str, subset_sequence_ids: list):
    """
    Given an input FastA file, subsets the file to the provided sequence IDs.

    :param input_fasta_filepath: Path to the input FastA file
    :param subset_sequence_ids: list of the names of sequences to keep. If any names in the list are not in the
                                input file, the function will not return anything for those names. The function will
                                raise an error if any duplicate sequence names are detected in the input FastA file.
    :return: generator of a SeqRecord object for the subset sequences
    """

    sequence_names = []
    for record in load_fasta_sequences(input_fasta_filepath):
        sequence_names.append(record.name)

        if record.name in subset_sequence_ids:
            yield record

    # Raise an error if there are duplicate sequence names
    if len(set(sequence_names)) < len(sequence_names):
        sequence_names_series = pd.Series(sequence_names)
        duplicates_names = set(sequence_names_series[sequence_names_series.duplicated() == True])

        error = RuntimeError(f'Duplicate sequence IDs were detected in the input FastA file "{input_fasta_filepath}": '
                             f'{", ".join(duplicates_names)}')
        logger.error(error)
        raise error


def generate_bed_file(contig_seqrecord: SeqIO.SeqRecord, bed_filepath: str, length_threshold: int = 100000):
    """
    Generates a BED file for a desired region around the ends of a contig. Writes the BED file to the bed_filepath.

    :param contig_seqrecord: SeqRecord of the contig
    :param bed_filepath: desired output filepath for the BED file
    :param length_threshold: length (bp) around the contig end to target in the BED file.
                             Half of this length will be returned around each end of the contig.
    """

    contig_name = contig_seqrecord.name
    contig_length = len(contig_seqrecord.seq)

    logger.debug('Making BED file with end proximity threshold of ' + str(length_threshold))

    if contig_length < length_threshold:
        logger.warning(f'Contig length ({contig_length}) is less than the supplied length threshold '
                       f'({length_threshold}), so the BED file will be for the whole contig.')

        contigs = [contig_name]
        starts = [0]
        stops = [contig_length]
    else:
        half_threshold = int(length_threshold / 2)
        contigs = [contig_name, contig_name]
        starts = [0, contig_length - half_threshold]
        stops = [half_threshold, contig_length]

    end_regions = pd.DataFrame({'contig': contigs, 'start': starts, 'stop': stops})

    end_regions.to_csv(bed_filepath, sep='\t', header=None, index=False)


def link_contig_ends(contig_record: SeqIO.SeqRecord, bam_filepath: str, length_outdir: str, length_threshold: int,
                     cli_tool_settings_dict: dict, verbose_logfile: str, threads: int = 1):
    """
    Attempt to stitch the ends of an input circular contig via assembling reads mapped within x bp of the contig ends.

    :param contig_record: SeqRecord of the circular contig
    :param bam_filepath: path to a BAM file with mapping information of long reads to the contig (it is OK if this file
                         also contains mappings to other contigs outside the contig of interest)
    :param length_outdir: output directory for the analysis; will be created by the function, although it can exist
                          before the function is run
    :param length_threshold: bp region around the contig ends to subset for the assembly (half of this distance will be
                             targeted from each contig end)
    :param cli_tool_settings_dict: dictionary containing the following CLI tool settings, as defined in main(), as keys:
                                   flye_read_mode, flye_read_error, circlator_min_id, circlator_min_length,
                                   circlator_ref_end, circlator_reassemble_end
    :param verbose_logfile: path to a log file where shell script logs will be saved
    :param threads: parallel processor threads to use for the analysis
    :return: exit status code for Flye
    """

    # Define folder and file names
    bed_filepath = os.path.join(length_outdir, 'ends.bed')
    ends_fastq_filepath = os.path.join(length_outdir, 'ends.fastq.gz')
    flye_length_outdir = os.path.join(length_outdir, 'assembly')
    merge_dir = os.path.join(length_outdir, 'merge')

    os.makedirs(length_outdir, exist_ok=True)

    # Get reads for the selected region around the contig ends
    generate_bed_file(contig_record, bed_filepath, length_threshold=length_threshold)
    subset_reads_from_bam(bam_filepath=bam_filepath, bed_filepath=bed_filepath,
                          subset_fastq_filepath=ends_fastq_filepath, log_filepath=verbose_logfile,
                          append_log=True, threads=threads)

    # Assemble the reads to get (hopefully) a joined contig end
    flye_exit_status = run_flye(fastq_filepath=ends_fastq_filepath, flye_outdir=flye_length_outdir,
                                flye_read_mode=cli_tool_settings_dict['flye_read_mode'],
                                flye_read_error=cli_tool_settings_dict['flye_read_error'],
                                log_filepath=verbose_logfile, append_log=True, threads=threads)

    if flye_exit_status == 0:
        # Write the original contig sequence to a file
        circular_contig_filepath = os.path.join(length_outdir, 'original.fasta')

        with open(circular_contig_filepath, 'w') as output_handle:
            SeqIO.write(contig_record, output_handle, 'fasta')

        # Stitch the joined contig end onto the original assembly
        # TODO - sometimes small contigs are already rotated far from original origin. Stitch point is
        #        hard to find. Does circlator report stitch point?
        run_circlator_merge(circular_contig_filepath=circular_contig_filepath,
                            patch_contig_filepath=os.path.join(flye_length_outdir, 'assembly.fasta'),
                            merge_outdir=merge_dir, circlator_min_id=cli_tool_settings_dict['circlator_min_id'],
                            circlator_min_length=cli_tool_settings_dict['circlator_min_length'],
                            circlator_ref_end=cli_tool_settings_dict['circlator_ref_end'],
                            circlator_reassemble_end=cli_tool_settings_dict['circlator_reassemble_end'],
                            log_filepath=verbose_logfile, append_log=True)
    else:
        logger.warning('Flye assembly FAILED')

    return flye_exit_status


def iterate_linking_contig_ends(contig_record: SeqIO.SeqRecord, bam_filepath: str, linking_outdir: str,
                                length_thresholds: list, cli_tool_settings_dict: dict, verbose_logfile: str,
                                threads: int = 1):
    """
    Iterate link_contig_ends to try to stitch the ends of a circular contig using multiple length thresholds.

    :param contig_record: SeqRecord of the circular contig
    :param bam_filepath: path to a BAM file with mapping information of long reads to the contig (it is OK if this file
                         also contains mappings to other contigs outside the contig of interest)
    :param linking_outdir: output directory for the analysis; will be created by the function, although it can exist
                           before the function is run
    :param length_thresholds: list of bp regions around the contig ends to attempt to subset for the assembly
    :param cli_tool_settings_dict: dictionary containing the following CLI tool settings, as defined in main(), as keys:
                                   flye_read_mode, flye_read_error, circlator_min_id, circlator_min_length,
                                   circlator_ref_end, circlator_reassemble_end
    :param verbose_logfile: path to a logfile where shell script logs will be saved
    :param threads: parallel processor threads to use for the analysis
    :return: boolean of whether end linkage was successful (True) or not (False)
    """

    os.makedirs(linking_outdir, exist_ok=True)

    # Where important log files will be copied
    log_dir_base = os.path.join(linking_outdir, 'logs')

    # Keep trying to link the contig ends until successful or until all length thresholds have been attempted
    assembly_attempts = 0
    linked_ends = False

    while linked_ends is False:
        # Exit the function if all length thresholds have been tried
        if assembly_attempts >= len(length_thresholds):
            break

        # Get the length threshold for this attempt
        length_threshold = length_thresholds[assembly_attempts]

        # Don't attempt to link the ends if the contig is too short
        if len(contig_record.seq) <= length_threshold:
            logger.info(f'Skipping length threshold of {length_threshold} because '
                        f'contig is shorter than this length ({len(contig_record.seq)} bp)')
            assembly_attempts = assembly_attempts + 1
            continue

        # Make sure that the circlator_min_length is shorter than the length threshold (otherwise merging is impossible)
        if cli_tool_settings_dict['circlator_min_length'] > length_threshold:
            revised_min_length = int(length_threshold * 0.9)
            cli_tool_settings_dict_revised = cli_tool_settings_dict.copy()
            cli_tool_settings_dict_revised['circlator_min_length'] = revised_min_length

            logger.warning(f'The minimum required length of alignment between the original and reassembled contig '
                           f'(specified by circlator_min_length; {cli_tool_settings_dict["circlator_min_length"]} bp) '
                           f'is longer than the length_threshold the original contig will be subset to '
                           f'({length_threshold} bp). This means that finding a matching merge will be impossible. To '
                           f'overcome this, the script will shorter the circlator_min_length for this iteration to 90% '
                           f'of the length threshold, i.e., {revised_min_length} bp.')

            cli_tool_settings_dict_used = cli_tool_settings_dict_revised

        else:
            cli_tool_settings_dict_used = cli_tool_settings_dict

        # Try to stitch the contig ends
        logger.info(f'Attempting reassembly with a length threshold of {length_threshold} bp')
        length_outdir = os.path.join(linking_outdir, 'tmp', f'L{length_threshold}')

        flye_exit_status = link_contig_ends(contig_record=contig_record, bam_filepath=bam_filepath,
                                            length_outdir=length_outdir, length_threshold=length_threshold,
                                            cli_tool_settings_dict=cli_tool_settings_dict_used,
                                            verbose_logfile=verbose_logfile, threads=threads)

        # Make the output logging directory (this is needed even if Flye fails so the pipeline can keep going)
        log_dir = os.path.join(log_dir_base, f'L{length_threshold}')
        os.makedirs(log_dir, exist_ok=True)

        # Stop early if Flye fails
        if flye_exit_status != 0:
            assembly_attempts = assembly_attempts + 1
            continue

        # If Flye was successful, copy important log files for that specific length threshold
        shutil.copy(os.path.join(length_outdir, 'assembly', 'assembly_info.txt'), log_dir)
        shutil.copy(os.path.join(length_outdir, 'merge', 'merge.circularise.log'), log_dir)
        shutil.copy(os.path.join(length_outdir, 'merge', 'merge.circularise_details.log'), log_dir)

        # See if stitching the contigs ends worked
        if check_circlator_success(os.path.join(length_outdir, 'merge', 'merge.circularise.log')) is True:
            logger.info('Successfully linked contig ends')

            # Rotate to midpoint so that the stitched points can be polished more effectively downstream
            rotate_sequences_wf(fasta_filepath=os.path.join(length_outdir, 'merge', 'merge.fasta'),
                                output_filepath=os.path.join(linking_outdir, 'stitched.fasta'), rotate_type='fraction',
                                rotate_value_single=0.5, max_sequences_in_file=1, append=False, strip_descriptions=True)

            # Save a copy of the final circlator merge logfile in the main log directory
            shutil.copy(os.path.join(length_outdir, 'merge', 'merge.circularise_details.log'),
                        os.path.join(log_dir_base, f'{contig_record.name}_circlator_final.log'))
            linked_ends = True

        assembly_attempts = assembly_attempts + 1

    shutil.rmtree(os.path.join(linking_outdir, 'tmp'))

    return linked_ends


def stitch_all_contigs(repair_paths: RepairPaths, length_thresholds: list, cli_tool_settings_dict: dict, threads: int):
    """
    Run the iterate_linking_contig_ends function on all contigs in an input FastA file, i.e., attempt to stitch the ends
    of all the contigs (assumed circular) in the file. Writes stitched contigs to end_repaired_contigs_filepath.

    :param repair_paths: RepairPaths object containing paths to output files used in the repair process.
    :param length_thresholds: list of bp regions around the contig ends to attempt to subset for the assembly
    :param cli_tool_settings_dict: dictionary containing the following CLI tool settings, as defined in main(), as keys:
                                   flye_read_mode, flye_read_error, circlator_min_id, circlator_min_length,
                                   circlator_ref_end, circlator_reassemble_end
    :param threads: parallel processor threads to use for the analysis
    :return: list of the names of any contigs that could not be stitched successfully (list length will be zero if all
             contigs stitched successfully)
    """

    # Initialize the repaired contigs FastA file (so it will overwrite an old file rather than just append later)
    with open(repair_paths.end_repaired_contigs_filepath, 'w') as output_handle:
        output_handle.write('')

    # Make the tmp directory for output files and the tmp directory where final log files (if any) will be moved
    os.makedirs(os.path.join(repair_paths.linking_outdir_base, 'log_summary'), exist_ok=True)

    failed_contig_names = []
    with open(repair_paths.circular_contig_tmp_fasta) as fasta_handle:
        for contig_record in SeqIO.parse(fasta_handle, 'fasta'):
            logger.info(f'Starting end repair on contig: {contig_record.name}')

            # This is the temp folder that will be used for this contig during stitching
            linking_outdir = os.path.join(repair_paths.linking_outdir_base, contig_record.name)
            end_linkage_complete = iterate_linking_contig_ends(contig_record, repair_paths.bam_filepath, linking_outdir,
                                                               length_thresholds, cli_tool_settings_dict,
                                                               repair_paths.verbose_logfile, threads)

            if end_linkage_complete is False:
                logger.warning(f'Contig {contig_record.name}: FAILED to linked contig ends')
                os.makedirs(os.path.join(repair_paths.linking_outdir_base, 'troubleshooting'), exist_ok=True)
                shutil.move(os.path.join(linking_outdir, 'logs'),
                            os.path.join(repair_paths.linking_outdir_base, 'troubleshooting', contig_record.name))

                failed_contig_names.append(contig_record.name)
            elif end_linkage_complete is True:
                # Append the successful contig onto the main file
                with open(os.path.join(linking_outdir, 'stitched.fasta')) as input_handle:
                    with open(repair_paths.end_repaired_contigs_filepath, 'a') as append_handle:
                        append_handle.write(input_handle.read())

                shutil.move(os.path.join(linking_outdir, 'logs', f'{contig_record.name}_circlator_final.log'),
                            os.path.join(repair_paths.linking_outdir_base, 'log_summary', f'{contig_record.name}.log'))
                shutil.rmtree(linking_outdir)
            else:
                error = ValueError(f'end_linkage_complete should be boolean True or False; is actually '
                                   f'{end_linkage_complete}')
                logger.error(error)
                raise error

    return failed_contig_names


def run_end_repair(long_read_filepath: str, assembly_info: AssemblyInfo, output_dir: str, length_thresholds: list,
                   keep_failed_contigs: bool, cli_tool_settings_dict: dict, threads: int, threads_mem_mb: int):
    """
    Runs the end repair workflow.

    :param long_read_filepath: path to the QC-passing Nanopore reads, FastQ format; gzipped is OK
    :param assembly_info: AssemblyInfo object with files paths containing info about the input assembly
    :param output_dir: path to the directory to save output files
    :param length_thresholds: list of bp regions around the contig ends to attempt to subset for the assembly
    :param keep_failed_contigs: boolean that defines whether to 1) continue the code even if some contigs cannot be end
                                 repaired (True) vs. to 2) exit with an error code if some contigs cannot be end
                                 repaired (False)
    :param cli_tool_settings_dict: dictionary containing the following CLI tool settings, as defined in main(), as keys:
                                   flye_read_mode, flye_read_error, circlator_min_id, circlator_min_length,
                                   circlator_ref_end, circlator_reassemble_end
    :param threads: parallel processor threads to use for the analysis
    :param threads_mem_mb: memory in MB per thread (to use for samtools); must be an integer
    """

    repair_paths = RepairPaths(output_dir)

    # TODO: for linear contigs, consider just getting all non-circular contigs regardless of whether they are in the
    #  assembly info file. This edit might require modification of the subset_sequences function and refactoring of code
    #  in several places, but I imagine it is closer to the behaviour that the user would expect.
    circular_contig_names, linear_contig_names = parse_assembly_info_file(assembly_info.assembly_info_filepath,
                                                                          assembly_info.assembly_info_type)

    # No need to run the pipeline if there are no circular contigs
    if len(circular_contig_names) == 0:
        repair_no_circular_input_contigs(assembly_info, repair_paths, linear_contig_names)

    # Subset circular contigs to use in the end repair workflow
    os.makedirs(repair_paths.linking_outdir_base, exist_ok=True)
    with open(repair_paths.circular_contig_tmp_fasta, 'w') as output_handle:
        for record in subset_sequences(assembly_info.assembly_fasta_filepath, circular_contig_names):
            SeqIO.write(record, output_handle, 'fasta')

    # Start the main workflow
    logger.info('Mapping reads to all contigs')
    map_long_reads(contig_filepath=assembly_info.assembly_fasta_filepath, long_read_filepath=long_read_filepath,
                   output_bam_filepath=repair_paths.bam_filepath, log_filepath=repair_paths.verbose_logfile,
                   append_log=False, threads=threads, threads_mem_mb=threads_mem_mb)

    failed_contig_names = stitch_all_contigs(repair_paths=repair_paths, length_thresholds=length_thresholds,
                                             cli_tool_settings_dict=cli_tool_settings_dict, threads=threads)
    os.makedirs(repair_paths.circlator_logs, exist_ok=True)
    shutil.move(os.path.join(repair_paths.linking_outdir_base, 'log_summary'), repair_paths.circlator_logs)

    # Check for contigs that could not be circularized
    if len(failed_contig_names) != 0:
        repair_handle_failed_contigs(assembly_info, repair_paths, failed_contig_names, keep_failed_contigs, output_dir)

    # Append linear contigs to the repaired contig file
    with open(repair_paths.end_repaired_contigs_filepath, 'a') as append_handle:
        for record in subset_sequences(assembly_info.assembly_fasta_filepath, linear_contig_names):
            SeqIO.write(record, append_handle, 'fasta')

    # Make a summary report
    contig_info = ContigInfo(circular_contig_names, failed_contig_names, linear_contig_names)
    write_repair_report(contig_info, repair_paths)

    # Clean up temp files
    os.remove(repair_paths.bam_filepath)
    os.remove(f'{repair_paths.bam_filepath}.bai')
    shutil.rmtree(repair_paths.linking_outdir_base)

    log_repair_complete(contig_info, repair_paths)


def log_repair_complete(contig_info, repair_paths):
    """
    Logs the completion of the repair process.

    :param contig_info: ContigInfo object containing lists of the names of repair-passed, repair-failed, and linear
                        contigs.
    :param repair_paths: RepairPaths object containing paths to output files used in the repair process.
    """
    logger.info('End repair finished!')
    logger.info('#### Final stats: ####')
    logger.info(f'Circular, repaired: {len(contig_info.repaired_contig_names)} contig(s)')
    logger.info(f'Circular, not repairable: {len(contig_info.failed_contig_names)} contig(s)')
    logger.info(f'Linear: {len(contig_info.linear_contig_names)} contig(s)')
    logger.info('######################')
    logger.info(f'Output contigs are saved at {repair_paths.end_repaired_contigs_filepath}. '
                f'A summary of repair work is saved at {repair_paths.end_repair_report_filepath}.')


def write_repair_report(contig_info, repair_paths):
    """
    This method writes a report of repair information to a TSV file.

    :param contig_info: ContigInfo object containing lists of the names of repair-passed, repair-failed, and linear
                        contigs.
    :param repair_paths: RepairPaths object containing paths to output files used in the repair process.
    """
    # Write info file of how contigs were repaired
    repair_info = pd.concat(
        [pd.DataFrame({'contig': contig_info.repaired_contig_names, 'circular': 'Y', 'repaired': 'Y'}),
         pd.DataFrame({'contig': contig_info.failed_contig_names, 'circular': 'Y', 'repaired': 'N'}),
         pd.DataFrame({'contig': contig_info.linear_contig_names, 'circular': 'N', 'repaired': 'N'})], axis=0)
    repair_info.to_csv(repair_paths.end_repair_report_filepath, sep='\t', index=False)


def repair_handle_failed_contigs(assembly_info: AssemblyInfo, repair_paths: RepairPaths, failed_contig_names,
                                 keep_failed_contigs, output_dir):
    """
    Repair and handle failed contigs in the assembly.

    :param assembly_info: AssemblyInfo object with files paths containing info about the input assembly.
    :param repair_paths: RepairPaths object containing paths to output files used in the repair process.
    :param failed_contig_names: A list of names of contigs that could not be circularized.
    :param keep_failed_contigs: A boolean indicating whether to keep the original (non-repaired) versions of
                                failed contigs in the final output file.
    :param output_dir: The directory to which the output files should be saved.
    """
    if os.path.isdir(os.path.join(repair_paths.linking_outdir_base, 'troubleshooting')):
        shutil.move(os.path.join(repair_paths.linking_outdir_base, 'troubleshooting'),
                    os.path.join(output_dir, 'troubleshooting'))
    if keep_failed_contigs is False:
        logger.error(f'{len(failed_contig_names)} contigs could not be circularized. A partial output file '
                     f'including successfully circularized contigs (and no linear contigs) is available at '
                     f'{repair_paths.end_repaired_contigs_filepath} for debugging. Exiting with error status. See temporary '
                     f'files and verbose logs for more details.')
        sys.exit(1)
    elif keep_failed_contigs is True:
        logger.warning(f'{len(failed_contig_names)} contigs could not be circularized. The original (non-repaired) '
                       f'versions of these contigs will be included in the final output file')
        logger.warning(f'Names of contigs that could not be circularized: {", ".join(failed_contig_names)}')

        # Get the non-repaired circular contigs and append them to the repaired contigs file
        with open(repair_paths.end_repaired_contigs_filepath, 'a') as append_handle:
            for record in subset_sequences(assembly_info.assembly_fasta_filepath, failed_contig_names):
                SeqIO.write(record, append_handle, 'fasta')
    else:
        error = ValueError(f'keep_failed_contigs should be a boolean True or False; you provided '
                           f'{keep_failed_contigs}')
        logger.error(error)
        raise error


def repair_no_circular_input_contigs(assembly: AssemblyInfo, repair_paths: RepairPaths, linear_contig_names: list):
    """
    Handle repair when no circular contigs are found.

    :param assembly: AssemblyInfo object with files paths containing info about the input assembly.
    :param repair_paths: RepairPaths object containing paths to output files used in the repair process.
    :param linear_contig_names: list of the names of non-circular contigs in the assembly info file
    """
    logger.info('No circular contigs. Will copy the input file, repaired_info.tsv, and finish early.')
    shutil.copyfile(assembly.assembly_fasta_filepath, repair_paths.end_repaired_contigs_filepath)

    repair_info = pd.DataFrame({'contig': linear_contig_names, 'circular': 'N', 'repaired': 'N'})
    repair_info.to_csv(repair_paths.end_repair_report_filepath, sep='\t', index=False)
    logger.info('Pipeline finished.')
    sys.exit(0)


def subparse_cli(subparsers, parent_parser: argparse.ArgumentParser = None):
    """
    Parses the CLI arguments and adds them as a subparser to an existing parser.

    :param subparsers: A special subparser action object created from an existing parser by the add_subparsers() method.
                       For example, parser = argparse.ArgumentParser(); subparsers = parser.add_subparsers().
    :param parent_parser: An optional ArgParse object with additional arguments (e.g., shared across all modules) to
                          add to this CLI parser. This can be a unique parser and does not need to be the parent of the
                          subparsers object. If None, then no parent will be added for the subparser.
    :return: An ArgumentParser object created by subparsers.add_parser()
    """

    description = 'Repairs ends of circular contigs from Flye.'

    # Initialize within the provided subparser
    subparser = subparsers.add_parser('repair', help=description, parents=[parent_parser] if parent_parser else [])

    # Add attribute to tell main() what sub-command was called.
    subparser.set_defaults(repair=True)

    required_settings = subparser.add_argument_group('Required')
    read_settings = subparser.add_argument_group('Input read options')
    merge_settings = subparser.add_argument_group('Merge options')
    workflow_settings = subparser.add_argument_group('Workflow options')

    required_settings.add_argument('-l', '--long_read_filepath', metavar='PATH', required=True, type=str,
                                   help='QC-passing Nanopore reads')
    required_settings.add_argument('-a', '--assembly_fasta_filepath', metavar='PATH', required=True, type=str,
                                   help='Contigs to be end-repaired')
    required_settings.add_argument('-i', '--assembly_info_filepath', metavar='PATH', required=True, type=str,
                                   help='assembly_info.txt file from Flye showing which assembled contigs are circular '
                                        'vs. linear (or a custom guide file; see -c below)')
    required_settings.add_argument('-o', '--output_dir', metavar='PATH', required=True, type=str,
                                   help='Output directory path')

    read_settings.add_argument('-f', '--flye_read_mode', metavar='MODE', required=False, default='nano-hq',
                               choices=['nano-hq', 'nano-raw'],
                               help='Type of long read reads provided by -l, to be used for reassembly by Flye. See '
                                    'details on these settings in the Flye documentation. (default: nano-hq)')
    read_settings.add_argument('-F', '--flye_read_error', metavar='FRACTION', required=False, default=0, type=float,
                               help='Expected error rate of input reads, expressed as proportion (e.g., 0.03). '
                                    'If "0", then have flye set the read error automatically (default: 0)')

    merge_settings.add_argument('-I', '--circlator_min_id', metavar='PERCENT', required=False, default=99, type=float,
                                help='Percent identity threshold for circlator merge (default: 99)')
    merge_settings.add_argument('-L', '--circlator_min_length', metavar='LENGTH', required=False, default=10000,
                                type=int, help='Minimum required overlap (bp) between original and merge contigs '
                                               '(default: 10000)')
    merge_settings.add_argument('-e', '--circlator_ref_end', metavar='LENGTH', required=False, default=100, type=int,
                                help='Minimum distance (bp) between end of original contig and nucmer hit '
                                     '(default: 100)')
    merge_settings.add_argument('-E', '--circlator_reassemble_end', metavar='LENGTH', required=False, default=100,
                                type=int, help='Minimum distance (bp) between end of merge contig and nucmer hit '
                                               '(default: 100)')

    workflow_settings.add_argument('-k', '--keep_going_with_failed_contigs', required=False, action='store_true',
                                   help='Set this flag to continue running this script even if some contigs '
                                        'cannot be circularized and end-repaired. For any non-repaired contigs, an '
                                        'exact copy of the original contig will be output.')
    workflow_settings.add_argument('-c', '--custom_assembly_info_file', required=False, action='store_true',
                                   help='Set this flag if you provide a custom tab-separated file to specify the '
                                        'circular vs. linear status of each contig in the assembly, rather than using '
                                        'the assembly_info.txt file output by Flye (default), as your argument to -i. '
                                        'The custom file must have the following format: no headers; tab-separated; '
                                        'first column is the contig names; second column is the status of the contigs, '
                                        'either "circular" or "linear". Any contig names not in this file will be '
                                        'dropped by the script!')
    workflow_settings.add_argument('-T', '--length_thresholds', metavar='LIST', required=False,
                                   default='100000,75000,50000,25000,5000,2500,1000', type=str,
                                   help='Comma-separated list of length thresholds for reassembly around the contig '
                                        'ends (bp) (default: 100000,75000,50000,25000,5000,2500,1000)')
    workflow_settings.add_argument('-t', '--threads', metavar='JOBS', required=False, default=1, type=int,
                                   help='Number of processors threads to use (default: 1)')
    workflow_settings.add_argument('-m', '--threads_mem', metavar='GB', required=False, default=1, type=float,
                                   help='Memory (GB) to use **per thread** for samtools sort (default: 1)')

    return subparser
