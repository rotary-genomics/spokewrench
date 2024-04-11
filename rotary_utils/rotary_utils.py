#!/usr/bin/env python

"""
Description: A command-line interface for the Rotary-utils toolkit.
Copyright: Jackson M. Tsuji and Lee H. Bergstrand 2024
"""

import logging
import argparse

import rotary_utils.repair as repair

# Set up the logger
logger = logging.getLogger(__name__)
formatter = logging.Formatter('[ %(asctime)s ]: %(levelname)s: %(funcName)s: %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(formatter)
logger.addHandler(stream_handler)


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

    # Select the sub-command to run.
    if hasattr(args, 'repair'):
        repair.main(args)
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

    # Declare sub-commands
    repair.parse_cli(subparsers)

    return parser


if __name__ == '__main__':
    main()
