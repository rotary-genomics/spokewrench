#!/usr/bin/env python
# rotate.py
# Rotates circular DNA elements
# Copyright Jackson M. Tsuji and Lee H. Bergstrand 2024

import argparse
import logging
import os
import sys
import math

from Bio import SeqIO

from rotary_utils.utils import set_write_mode, load_fasta_sequences

logger = logging.getLogger(__name__)


# TODO: this is just a placeholder function for now to show some of the features I'd like to add to the CLI later on.
def main(args):
    """
    Starts the rotate utility.

    :param args: args parsed by parse_cli()
    """

    # Parse some command line inputs further
    # TODO: improve error handling if a string is provided instead of a real length
    # TODO: confirm what the default value of args.contig_names will be - might need to specify
    if args.contig_names is not None:
        contig_names = [int(x) for x in args.contig_names.split(',')]

    # Startup messages
    logger.info('Running ' + os.path.basename(sys.argv[0]))
    logger.debug('### SETTINGS ###')
    logger.debug(f'Input FastA: {args.input_fasta}')
    logger.debug(f'Output FastA: {args.output_fasta}')
    # TODO: finish the rest of the startup messages

    logger.debug(f'Threads: {args.threads}')
    logger.debug(f'Verbose logging: {args.verbose}')
    logger.debug('################')

    # TODO: add rotate code

    logger.info(os.path.basename(sys.argv[0]) + ': done.')


def parse_gene_positions():
    """
    Parses the positions of genes along a contig.

    :return: table of gene positions
    """
    # TODO: this function is a stub.

    pass


def avoid_gene_collision():
    """
    Checks if the desired rotation point will collide with a gene and provides an alternative safe point if needed.

    :return: The next-nearest position without a collision, with a buffer or on a particular side if requested
    """
    # TODO: this function is a stub.

    pass


def rotate_sequence(sequence_record: SeqIO.SeqRecord, rotate_position: int, strip_description: bool = True):
    """
    Rotates an input (circular) sequence to a specified position.

    :param sequence_record: BioPython SeqRecord object containing the sequence to rotate
    :param rotate_position: new start position for the contig after rotation (if 0, contig is not rotated)
    :param strip_description: boolean of whether to trim off the read description in the output sequences
    :return: BioPython SeqRecord object of the rotated FastA file
    """

    sequence_length = len(sequence_record.seq)
    if rotate_position > sequence_length:
        logger.error(f'Desired rotate position ({rotate_position}) is larger than the sequence length '
                     f'({sequence_length} bp).')
        raise RuntimeError

    elif rotate_position == sequence_length:
        logger.warning(f'Desired rotate position ({rotate_position}) is equal to the sequence length '
                       f'({sequence_length} bp), so the contig will effectively not be rotated.')

    elif rotate_position < 0:
        logger.debug(f'Desired rotate position ({rotate_position}) is less than 0.')
        raise RuntimeError

    logger.debug(f'Rotating sequence to start at position: {rotate_position}')
    contig_sequence_rotated_front = sequence_record.seq[rotate_position:sequence_length]
    contig_sequence_rotated_back = sequence_record.seq[0:rotate_position]
    contig_sequence_rotated = contig_sequence_rotated_front + contig_sequence_rotated_back

    # Update SeqRecord
    sequence_record.seq = contig_sequence_rotated

    if strip_description:
        sequence_record.description = sequence_record.name

    return sequence_record


def rotate_sequence_to_fraction(sequence_record: SeqIO.SeqRecord, rotate_fraction: float,
                                strip_description: bool = True):
    """
    Rotates an input (circular) sequence to a specified fraction of the total length.

    :param sequence_record: BioPython SeqRecord object containing the sequence to rotate
    :param rotate_fraction: fraction of the contig length to rotate the contig to (rounded down to nearest bp)
    :param strip_description: boolean of whether to trim off the read description in the output sequences
    :return: BioPython SeqRecord object of the rotated FastA file
    """

    sequence_length = len(sequence_record.seq)

    if rotate_fraction < 0:
        logger.error(f'Requested rotation fraction ({rotate_fraction}) is less than 0.')
        raise RuntimeError

    else:
        rotate_position = math.floor(sequence_length * rotate_fraction)

    sequence_record = rotate_sequence(sequence_record, rotate_position=rotate_position,
                                      strip_description=strip_description)

    return sequence_record


def rotate_sequence_to_midpoint(sequence_record: SeqIO.SeqRecord, strip_description: bool = True):
    """
    Rotates an input (circular) sequence to the approximate midpoint.

    :param sequence_record: BioPython SeqRecord object containing the sequence to rotate
    :param strip_description: boolean of whether to trim off the read description in the output sequences
    :return: BioPython SeqRecord object of the rotated FastA file
    """

    sequence_record = rotate_sequence_to_fraction(sequence_record, rotate_fraction=0.5,
                                                  strip_description=strip_description)
    return sequence_record


def rotate_sequences_wf(fasta_filepath: str, output_filepath: str, rotate_positions: dict, append: bool = False,
                        strip_descriptions: bool = True):
    """
    Rotates all (circular) sequences in an input FastA file to a specified position.

    :param fasta_filepath: Path to the FastA file (unzipped) to load. Assumes all entries are circular.
    :param output_filepath: path where the FastA file containing the output rotated contig should be saved
    :param rotate_positions: dict of contig names and new start positions for the contigs after rotation
                             (if 0, contig is not rotated)
    :param append: whether to append the output FastA onto an existing file (True) or overwrite (False)
    :param strip_descriptions: boolean of whether to trim off the read description in the output sequences
    :return: writes file to disk
    """

    sequence_records = load_fasta_sequences(fasta_filepath)

    sequence_records_rotated = []
    for record in sequence_records:
        sequence_record_rotated = rotate_sequence(record, rotate_position=rotate_positions[record.name],
                                                  strip_description=strip_descriptions)
        sequence_records_rotated.append(sequence_record_rotated)

    write_mode = set_write_mode(append)
    with open(output_filepath, write_mode) as output_handle:
        for record in sequence_records_rotated:
            SeqIO.write(record, output_handle, 'fasta')


def rotate_sequences_to_midpoint_wf(fasta_filepath: str, output_filepath: str, append: bool = False,
                                    strip_descriptions: bool = True):
    """
    Rotates all (circular) sequences in an input FastA file to their midpoint.

    :param fasta_filepath: Path to the FastA file (unzipped) to load. Assumes all entries are circular.
    :param output_filepath: path where the FastA file containing the output rotated contig should be saved
    :param append: whether to append the output FastA onto an existing file (True) or overwrite (False)
    :param strip_descriptions: boolean of whether to trim off the read description in the output sequences
    :return: writes file to disk
    """

    sequence_records = load_fasta_sequences(fasta_filepath)

    sequence_records_rotated = []
    for record in sequence_records:
        sequence_record_rotated = rotate_sequence_to_midpoint(record, strip_description=strip_descriptions)
        sequence_records_rotated.append(sequence_record_rotated)

    write_mode = set_write_mode(append)
    with open(output_filepath, write_mode) as output_handle:
        for record in sequence_records_rotated:
            SeqIO.write(record, output_handle, 'fasta')


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

    description = """Rotates circular contigs."""

    # Initialize within the provided subparser
    subparser = subparsers.add_parser('rotate', help=description, parents=[parent_parser] if parent_parser else [])

    # Add attribute to tell main() what sub-command was called.
    subparser.set_defaults(rotate=True)

    required_settings = subparser.add_argument_group('Required')
    rotate_settings = subparser.add_argument_group('Rotate options')
    workflow_settings = subparser.add_argument_group('Workflow options')

    required_settings.add_argument('-i', '--input_fasta', required=True, type=str,
                                   help='Input contig fasta file')
    required_settings.add_argument('-o', '--output_fasta', required=True, type=str,
                                   help='Output contig fasta file. (Note that all input contigs will always be '
                                        'written to output. If you select no rotate options, the input contigs will '
                                        'be output with no sequence changes.)')

    rotate_settings.add_argument('-m', '--midpoint', required=False, action='store_true',
                                 help='Rotate all contigs to their midpoint (overrides -p, -P, -t, and -T)')
    rotate_settings.add_argument('-p', '--rotate_position', required=False, type=int,
                                 help='Base pair position to rotate contigs to (to specify a unique rotate '
                                      'position for each contig, see -t). Incompatible with -P.')
    rotate_settings.add_argument('-t', '--rotate_position_table', required=False, type=str,
                                 help='Path to a tab-separated file that contains the exact rotate position of '
                                      'each contig. Columns (with headers) should be "contig_id" and "rotate_position". '
                                      'Overrides -p and -P. Only contigs specified in the table will be rotated, although '
                                      'all contigs in the input file will be written to the output file. Incompatible with '
                                      '-T.')
    rotate_settings.add_argument('-P', '--rotate_proportion', required=False, type=float,
                                 help='Fractional position to rotate contigs to (e.g., 0.3 to rotate 30% of total length). '
                                      'To specify a unique fractional position to rotate for each contig, see -T). '
                                      'Incompatible with -p.')
    rotate_settings.add_argument('-T', '--rotate_proportion_table', required=False, type=str,
                                 help='Path to a tab-separated file that contains the exact fractional positions to rotate '
                                      'each contig to. Columns (with headers) should be "contig_id" and "rotate_proportion". '
                                      'Overrides -p and -P. Only contigs specified in the table will be rotated, although '
                                      'all contigs in the input file will be written to the output file. Incompatible '
                                      'with -t.')
    rotate_settings.add_argument('-c', '--contig_names', required=False, type=str,
                                 help='Contigs to be rotated (comma-separated list of IDs). '
                                      'Rotate operations will only be applied to these contigs, '
                                      'although all contigs in the input file will be written to output. '
                                      'If used with -t or -T, the contigs listed in the table will be further subset '
                                      'to those that are also in the provided list of contig names.')
    rotate_settings.add_argument('-r', '--output_report', required=False, type=str,
                                 help='Path to write an optional report file that shows how contigs were rotated')
    rotate_settings.add_argument('-b', '--base_start_position', required=False, default=1, type=int,
                                 choices=[0, 1],
                                 help='Whether the first base position on a contig is considered as 1 (default) or 0.')

    workflow_settings.add_argument('-t', '--threads', required=False, default=1, type=int,
                                   help='Number of processors threads to use (default: 1)')
    workflow_settings.add_argument('-v', '--verbose', required=False, action='store_true',
                                   help='Enable verbose logging to screen')

    return subparser
