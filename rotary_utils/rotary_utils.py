#!/usr/bin/env python

"""
Description: A command-line interface for the Rotary-utils toolkit.
Copyright: Jackson M. Tsuji and Lee H. Bergstrand 2024
"""

import os
import logging
import argparse

from rotary_utils.utils import check_output_file, set_up_output_directory
from rotary_utils.rotate import main as rotate_main
from rotary_utils.rotate import subparse_cli as subparse_rotate_cli
from rotary_utils.repair import main as repair_main
from rotary_utils.repair import subparse_cli as subparse_repair_cli


def main():
    """
    Collects input arguments and selects a command to perform.
    """

    """
    When installing through pip with pyprojects.toml a new python script is generated
    that calls main() out of rotary_utils.py. Run parser inside main() so it can be called 
    externally as a function.
    """
    parser = parse_cli()
    args = parser.parse_args()

    # Initialize the root logger with a stream handler
    logger = logging.getLogger()
    formatter = logging.Formatter('[ %(asctime)s ]: %(levelname)s: %(filename)s: %(funcName)s: %(message)s',
                                  datefmt='%Y-%m-%d %H:%M:%S')
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    if (hasattr(args, 'verbose')) and (args.verbose is True):
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)

    # File logger setup
    if (hasattr(args, 'logfile')) and (args.logfile is not None):
        check_output_file(output_filepath=args.logfile, overwrite=args.overwrite)

        # Note: The logging level set here is the minimum level that the logger can write to. The actual level is
        #       defined for the root logger above in the args.verbose conditional.
        custom_path_file_handler = logging.FileHandler(filename=args.logfile, mode='w')
        custom_path_file_handler.setFormatter(formatter)
        custom_path_file_handler.setLevel(logging.DEBUG)
        logger.addHandler(custom_path_file_handler)

    if hasattr(args, 'output_dir'):
        set_up_output_directory(output_directory_filepath=args.output_dir, overwrite=args.overwrite)

        # Start log file in the output dir
        # Note: The logging level set here is the minimum level that the logger can write to. The actual level is
        #       defined for the root logger above in the args.verbose conditional.
        output_dir_file_handler = logging.FileHandler(filename=os.path.join(args.output_dir, 'log.txt'), mode='w')
        output_dir_file_handler.setFormatter(formatter)
        output_dir_file_handler.setLevel(logging.DEBUG)
        logger.addHandler(output_dir_file_handler)

    # Select the sub-command to run.
    if hasattr(args, 'rotate'):
        rotate_main(args)
    if hasattr(args, 'repair'):
        repair_main(args)
    else:
        parser.print_help()


def parse_cli():
    """
    Parses the CLI arguments.
    :return: An argparse parser object.
    """
    cli_title = """Utilities for manipulating genome or metagenome data with circular DNA elements. These
                   utilities integrate with the Rotary pipeline."""
    parser = argparse.ArgumentParser(description=cli_title)
    subparsers = parser.add_subparsers(help='Available Sub-commands')

    # Common arguments across all modules
    parent_parser = argparse.ArgumentParser(add_help=False)
    basic_config = parent_parser.add_argument_group('Basic config settings')
    basic_config.add_argument('-lf', '--logfile', metavar='PATH', required=False, default=None, type=str,
                              help='Log filepath (default: None)')
    basic_config.add_argument('-O', '--overwrite', required=False, action='store_true',
                              help='Overwrite existing files/directories. By setting this flag, you risk erasing old '
                                   'data.')
    basic_config.add_argument('-v', '--verbose', required=False, action='store_true',
                              help='Enable verbose logging')

    # Declare sub-commands
    subparse_rotate_cli(subparsers, parent_parser)
    subparse_repair_cli(subparsers, parent_parser)

    return parser


if __name__ == '__main__':
    main()
