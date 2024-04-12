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
import pandas as pd

from rotary_utils.utils import check_output_file, set_write_mode, load_fasta_sequences

logger = logging.getLogger(__name__)


def main(args):
    """
    Starts the rotate utility.

    :param args: args parsed by parse_cli()
    """

    # Startup messages
    logger.info('Running ' + os.path.basename(sys.argv[0]))
    logger.debug('### SETTINGS ###')
    logger.debug(f'Input FastA: {args.input_fasta}')
    logger.debug(f'Output FastA: {args.output_fasta}')
    logger.debug(f'Midpoint rotation: {args.midpoint}')
    logger.debug(f'Constant rotate position: {args.rotate_position}')
    logger.debug(f'Constant rotate fraction: {args.rotate_fraction}')
    logger.debug(f'Rotate position table: {args.rotate_position_table}')
    logger.debug(f'Rotate fraction table: {args.rotate_fraction_table}')
    logger.debug(f'Sequence names to rotate: {args.sequence_names}')
    logger.debug(f'Rotation report: {args.output_report}')
    logger.debug(f'Strip FastA header descriptions: {args.strip_descriptions}')
    logger.debug(f'Overwrite output files if needed?: {args.overwrite}')
    logger.debug(f'Verbose logging: {args.verbose}')
    logger.debug('################')

    # Check that the user has provided OK instructions
    check_rotate_args(args)

    # TODO: improve error handling if a string is provided instead of a real length
    if args.sequence_names is not None:
        sequence_names = [int(x) for x in args.sequence_names.split(',')]
        logger.debug(f'User provided a list of {len(args.sequence_names)} sequence IDs to search for rotation.')
    else:
        sequence_names = None

    # Perform the rotate operation based on which of the rotate arguments was specified by the user
    if args.midpoint is True:
        report = rotate_sequences_wf(args.input_fasta, args.output_fasta, rotate_type='fraction',
                                     rotate_value_single=0.5, sequence_names=sequence_names,
                                     strip_descriptions=args.strip_descriptions)

    elif args.rotate_position is not None:
        report = rotate_sequences_wf(args.input_fasta, args.output_fasta, rotate_type='position',
                                     rotate_value_single=args.rotate_position, sequence_names=sequence_names,
                                     strip_descriptions=args.strip_descriptions)

    elif args.rotate_fraction is not None:
        report = rotate_sequences_wf(args.input_fasta, args.output_fasta, rotate_type='fraction',
                                     rotate_value_single=args.rotate_fraction, sequence_names=sequence_names,
                                     strip_descriptions=args.strip_descriptions)

    elif args.rotate_position_table is not None:
        rotate_values_dict = parse_rotate_value_table(args.rotate_position_table, rotate_type='position')

        report = rotate_sequences_wf(args.input_fasta, args.output_fasta, rotate_type='position',
                                     rotate_values=rotate_values_dict, sequence_names=sequence_names,
                                     strip_descriptions=args.strip_descriptions)

    elif args.rotate_fraction_table is not None:
        rotate_values_dict = parse_rotate_value_table(args.rotate_position_table, rotate_type='fraction')

        report = rotate_sequences_wf(args.input_fasta, args.output_fasta, rotate_type='fraction',
                                     rotate_values=rotate_values_dict, sequence_names=sequence_names,
                                     strip_descriptions=args.strip_descriptions)

    else:
        raise RuntimeError('Could not find any rotate commands to execute.')

    # Save the rotation report if desired
    if args.output_report is not None:
        report.to_csv(args.output_report, sep='\t', index=False)

    logger.info(os.path.basename(sys.argv[0]) + ': done.')


def check_rotate_args(args):
    """
    Checks that the arguments supplied to the rotate module are non-conflicting. Errors if a conflict is found.

    :param args: Arguments parsed by subparse_cli()
    """

    """
    1. As described in the README, only one of the following arguments is allowed to be set for a run:
    args.midpoint
    args.rotate_position
    args.rotate_fraction
    args.rotate_position_table
    args.rotate_fraction_table
    
    2. If overwrite is False, then no output files can already exist:
    args.output_fasta
    args.output_report
    """

    # 1. Check rotate value arguments
    rotate_value_arguments = [args.midpoint, args.rotate_position, args.rotate_fraction, args.rotate_position_table,
                              args.rotate_fraction_table]

    defined_arguments = 0
    for argument in rotate_value_arguments:
        if argument is not None:
            defined_arguments = defined_arguments + 1

    if defined_arguments > 1:
        raise RuntimeError('More than one of -m, -p, -P, -t, and -T were specified in the command line, but only one '
                           'of these can be set for a rotate run.')
    elif defined_arguments == 0:
        raise RuntimeError('None of -m, -p, -P, -t, and -T were specified in the command line, but you must select one '
                           'of these to perform a rotate run.')
    elif defined_arguments != 0:
        raise RuntimeError('Ran into an issue processing the -m, -p, -P, -t, and -T arguments in the CLI.')

    # 2. Check output files
    check_output_file(output_filepath=args.output_fasta, overwrite=args.overwrite)

    if args.output_report is not None:
        check_output_file(output_filepath=args.output_report, overwrite=args.overwrite)


def parse_rotate_value_table(rotate_table_filepath: str, rotate_type: str):
    """
    Parses a tab-separated table as defined in the CLI containing rotation information.

    :param rotate_table_filepath: path to the tab-separated table with sequence rotation info
    :param rotate_type: either 'position' or 'fraction'. Runs a check to make sure column names match the specified
                        rotation type.
    :return: dictionary of sequence names (keys) and rotation values (values)
    """

    if rotate_type == 'position':
        expected_column_names = ['sequence_id', 'rotate_position']
    elif rotate_type == 'fraction':
        expected_column_names = ['sequence_id', 'rotate_fraction']
    else:
        raise ValueError(f'Expected "position" or "fraction", but got "{rotate_type}".')

    rotate_table = pd.read_csv(rotate_table_filepath, sep='\t')

    if expected_column_names not in rotate_table.columns:
        raise RuntimeError(f'Because rotate_type was set as "{rotate_type}", expected the column names '
                           f'"{",".join(expected_column_names)}", but could not find these in the loaded columns.'
                           f'Instead, found "{",".join(rotate_table.columns)}"')

    rotate_values_dict = dict(zip(rotate_table[expected_column_names[0]], rotate_table[expected_column_names[1]]))

    return rotate_values_dict


def parse_gene_positions():
    """
    Parses the positions of genes along a sequence.

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


def rotate_sequence_to_position(sequence_record: SeqIO.SeqRecord, rotate_position: int, strip_description: bool = True):
    """
    Rotates an input (circular) sequence to a specified position.

    :param sequence_record: BioPython SeqRecord object containing the sequence to rotate
    :param rotate_position: new start position for the sequence after rotation (if 0, sequence is not rotated)
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
                       f'({sequence_length} bp), so the sequence will effectively not be rotated.')

    elif rotate_position < 0:
        logger.debug(f'Desired rotate position ({rotate_position}) is less than 0.')
        raise RuntimeError

    logger.debug(f'Rotating sequence to start at position: {rotate_position}')
    sequence_rotated_front = sequence_record.seq[rotate_position:sequence_length]
    sequence_rotated_back = sequence_record.seq[0:rotate_position]
    sequence_rotated = sequence_rotated_front + sequence_rotated_back

    # Update SeqRecord
    sequence_record.seq = sequence_rotated

    if strip_description:
        sequence_record.description = sequence_record.name

    return sequence_record


def rotate_sequence_to_fraction(sequence_record: SeqIO.SeqRecord, rotate_fraction: float,
                                strip_description: bool = True):
    """
    Rotates an input (circular) sequence to a specified fraction of the total length.

    :param sequence_record: BioPython SeqRecord object containing the sequence to rotate
    :param rotate_fraction: fraction of the sequence length to rotate the sequence to (rounded down to nearest bp). Set
                            to 0.5 to rotate to midpoint.
    :param strip_description: boolean of whether to trim off the read description in the output sequences
    :return: BioPython SeqRecord object of the rotated FastA file
    """

    sequence_length = len(sequence_record.seq)

    if rotate_fraction < 0:
        logger.error(f'Requested rotation fraction ({rotate_fraction}) is less than 0.')
        raise RuntimeError

    else:
        rotate_position = math.floor(sequence_length * rotate_fraction)
        logger.debug(f'Rotating contig to fraction {rotate_fraction} = position {rotate_position} bp')

    sequence_record = rotate_sequence_to_position(sequence_record, rotate_position=rotate_position,
                                                  strip_description=strip_description)

    return sequence_record


def rotate_sequences_wf(fasta_filepath: str, output_filepath: str, rotate_type: str, rotate_values: dict = None,
                        rotate_value_single: float = None, sequence_names: list = None, max_sequences_in_file: int = 0,
                        append: bool = False, strip_descriptions: bool = True):
    """
    Rotates all (circular) sequences in an input FastA file to specified positions or fractions.

    :param fasta_filepath: Path to the FastA file (unzipped) to load
    :param output_filepath: path where the FastA file containing the output rotated sequence should be saved
    :param rotate_type: whether the rotation is a 'position' based or a 'fraction' based rotation
    :param rotate_values: dict of sequence names (keys) and the desired rotation position or fraction for each sequence
                          (values). If rotate_type if 'position', expects the values to be the new start position for
                          each sequence after rotation. If rotate_type if 'fraction', expects the values to be the
                          fraction of total sequence length to rotate. If the value is set to 0 for a sequence, then
                          that sequence is not rotated.
    :param rotate_value_single: a single value (position or fraction, depending on rotate_type) to rotate all sequences
                                to. This is useful if you want to rotate all sequences to their midpoint, for example,
                                via rotate_value_single=0.5 and rotate_type='fraction'. As a warning, if doing
                                positional rotation, the largest you can rotate to is the length of the shortest
                                sequence. This param cannot be set at the same time as rotate_values.
    :param sequence_names: list of sequence header IDs. Rotate operations will only be performed on these names.
    :param max_sequences_in_file: maximum number of allowable sequences in the input file (default: 0 = unlimited)
    :param append: whether to append the output FastA onto an existing file (True) or overwrite (False)
    :param strip_descriptions: boolean of whether to trim off the read description in the output sequences
    :return: tabular report of rotation stats for each sequence
    """

    subset_rotation = True if sequence_names is not None else False

    # Parse what kind of rotation values were supplied
    if (rotate_values is not None) and (rotate_value_single is not None):
        raise RuntimeError('Cannot set both rotate_values and rotate_value_single.')
    elif (rotate_values is None) and (rotate_value_single is None):
        raise RuntimeError('Must set one of rotate_values and rotate_value_single.')

    sequence_ids = []
    sequence_lengths = []
    final_rotate_values = []
    total_sequences_in_file = 0

    # Rotate each sequence by the specified amount
    write_mode = set_write_mode(append)
    with open(output_filepath, write_mode) as output_handle:
        for record in load_fasta_sequences(fasta_filepath):

            # If rotating specific sequences is turned off, then add the specified sequence ID to the rotation list so
            # that it passes the sequence name check
            if subset_rotation is False:
                sequence_names = [record.name]

            if record.name in sequence_names:

                # Get the rotation value for the sequence record
                if rotate_values is not None:
                    try:
                        rotate_value = rotate_values[record.name]
                    except KeyError:
                        logger.debug(f'Sequence name {record.name} was not supplied in rotate_values. Will not rotate.')
                        rotate_value = 0

                elif rotate_value_single is not None:
                    rotate_value = rotate_value_single

                else:
                    raise RuntimeError('Could not find rotate_values or rotate_value_single.')

                # Rotate the sequence
                if rotate_type == 'position':
                    record_rotated = rotate_sequence_to_position(record, rotate_position=rotate_value,
                                                                 strip_description=strip_descriptions)
                elif rotate_type == 'fraction':
                    record_rotated = rotate_sequence_to_fraction(record, rotate_fraction=rotate_value,
                                                                 strip_description=strip_descriptions)
                else:
                    raise ValueError(f'Expected "position" or "fraction", but not {rotate_type}')

            else:
                # Just pass the record through without rotation
                record_rotated = rotate_sequence_to_position(record, rotate_position=0,
                                                             strip_description=strip_descriptions)

            # Get values for output report
            sequence_ids.append(record.name)
            sequence_lengths.append(len(record.seq))
            final_rotate_values.append(rotate_value)

            total_sequences_in_file = total_sequences_in_file + 1

            if (max_sequences_in_file != 0) & (max_sequences_in_file < total_sequences_in_file):
                raise RuntimeError(f'A maximum of {max_sequences_in_file} sequences is allowed, but the file contains '
                                   f'more sequences than this.')

            SeqIO.write(record_rotated, output_handle, 'fasta')

    output_report = pd.DataFrame({'sequence_id': sequence_ids, 'sequence_length_bp': sequence_lengths,
                                  'rotate_type': rotate_type, 'rotate_value': final_rotate_values})

    return output_report


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

    description = """Rotates circular sequence."""

    # Initialize within the provided subparser
    subparser = subparsers.add_parser('rotate', help=description, parents=[parent_parser] if parent_parser else [])

    # Add attribute to tell main() what sub-command was called.
    subparser.set_defaults(rotate=True)

    required_settings = subparser.add_argument_group('Required')
    rotate_settings = subparser.add_argument_group('Rotate options')
    workflow_settings = subparser.add_argument_group('Workflow options')

    required_settings.add_argument('-i', '--input_fasta', metavar='PATH', required=True, type=str,
                                   help='Input fasta file')
    required_settings.add_argument('-o', '--output_fasta', metavar='PATH', required=True, type=str,
                                   help='Output fasta file. Note that all input sequences will always be '
                                        'written to output, even if only a subset of sequences are selected for '
                                        'rotate operations.')

    rotate_settings.add_argument('-m', '--midpoint', required=False, action='store_true',
                                 help='Rotate all sequences to their midpoint. (Incompatible with -p, -P, -t, and -T.)')
    rotate_settings.add_argument('-p', '--rotate_position', metavar='INT', required=False, type=int,
                                 help='Base pair position to rotate all sequences to (to specify a unique rotate '
                                      'position for each sequence, see -t). Incompatible with -m, -P, -t, and -T.')
    rotate_settings.add_argument('-P', '--rotate_fraction', metavar='FRACTION', required=False, type=float,
                                 help='Fractional position to rotate sequences to (e.g., 0.3 to rotate 30% of total '
                                      'length). To specify a unique fractional position to rotate for each sequence, '
                                      'see -T). Incompatible with -m, -p, -t, and -T.')
    rotate_settings.add_argument('-t', '--rotate_position_table', metavar='PATH', required=False, type=str,
                                 help='Path to a tab-separated file that contains the exact rotate position of '
                                      'each sequence. Columns (with headers) should be "sequence_id" and '
                                      '"rotate_position". Only sequences specified in the table will be rotated, '
                                      'although all sequences in the input file will be written to the output file. '
                                      'Incompatible with -m, -p, -P, and -T.')
    rotate_settings.add_argument('-T', '--rotate_fraction_table', metavar='PATH', required=False, type=str,
                                 help='Path to a tab-separated file that contains the exact fractional positions to '
                                      'rotate each sequence to. Columns (with headers) should be "sequence_id" and '
                                      '"rotate_fraction". Only sequences specified in the table will be rotated, '
                                      'although all sequences in the input file will be written to the output file. '
                                      'Incompatible with -m, -p, -P, and -t.')
    rotate_settings.add_argument('-n', '--sequence_names', metavar='LIST', required=False, type=str, default=None,
                                 help='Sequences to be rotated (comma-separated list of IDs). Rotate operations will '
                                      'only be applied to these sequences, although all sequences in the input file '
                                      'will be written to output. If used with -t or -T, the sequences listed in the '
                                      'table will be further subset to those that are also in the provided list of '
                                      'sequence names.')
    rotate_settings.add_argument('-r', '--output_report', metavar='PATH', required=False, type=str,
                                 help='Path to write an optional tab-separated report file that shows how sequences '
                                      'were rotated.')

    workflow_settings.add_argument('-s', '--strip_descriptions', required=False, action='store_true',
                                   help='Strip descriptions off FastA headers (i.e., any text after the first '
                                        'whitespace in each header)')

    return subparser
