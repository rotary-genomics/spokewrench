# spokewrench
Circularization-related utilities for the rotary pipeline

## Overview
*spokewrench*, currently under development, is a suite of utilities for manipulating circular DNA elements.
Once completed, these utilities will be able to serve as a modern replacement for
[circlator](https://github.com/sanger-pathogens/circlator), which is now in a frozen development state.
The utilities in *spokewrench* are built into the [*rotary* project](https://github.com/rotary-genomics/rotary),
where they are used with "best practices" to ensure proper assembly of circular sequences. These utilities can also
be used in a standalone fashion or integrated into custom genome assembly workflows to ensure circular DNA/RNA elements
are assembled accurately, similarly to how circlator is currently used.

## Requirements
- OS: Runs on Linux (tested on Ubuntu 20.04 and Ubuntu 22.04) and macOS (tested on Sonoma 14)
- Software: requires miniconda or manual installation using the dependencies shown in spokewrench/environment.yml
- Resources: should run on a modern laptop with >=8 GB RAM and >=4 CPU threads, in most cases

## Installation
```bash
git clone https://github.com/rotary-genomics/spokewrench.git

conda env create -n spokewrench --file=spokewrench/environment.yml
# Note: On Macs with Apple Silicon, you might need to add the flag '--platform osx-64' to force x64 install.
#       The installed packages will then be converted to ARM architecture using Rosetta on first use.

conda activate spokewrench

cd spokewrench

pip install --editable .

# See available commands
spokewrench -h
```

## Modules
### `repair`
Run on the outputs of Flye to repair the
[short gap or overlap region](https://github.com/fenderglass/Flye/issues/315#issuecomment-720679812) that can occur at
the ends of circular contigs produced by this assembly tool.

Help menu (`spokewrench repair -h`):
```commandline
usage: spokewrench repair [-h] [-lf PATH] [-O] [-v] -l PATH -a PATH -i PATH -o PATH [-f MODE] [-F FRACTION] [-I PERCENT] [-L LENGTH] [-e LENGTH] [-E LENGTH] [-k] [-c]
                          [-T LIST] [-t JOBS] [-m GB]

optional arguments:
  -h, --help            show this help message and exit

Basic config settings:
  -lf PATH, --logfile PATH
                        Log filepath (default: None)
  -O, --overwrite       Overwrite existing files/directories. By setting this flag, you risk erasing old data.
  -v, --verbose         Enable verbose logging

Required:
  -l PATH, --long_read_filepath PATH
                        QC-passing Nanopore reads
  -a PATH, --assembly_fasta_filepath PATH
                        Contigs to be end-repaired
  -i PATH, --assembly_info_filepath PATH
                        assembly_info.txt file from Flye showing which assembled contigs are circular vs. linear (or a custom guide file; see -c below)
  -o PATH, --output_dir PATH
                        Output directory path

Input read options:
  -f MODE, --flye_read_mode MODE
                        Type of long read reads provided by -l, to be used for reassembly by Flye. See details on these settings in the Flye documentation. (default: nano-hq)
  -F FRACTION, --flye_read_error FRACTION
                        Expected error rate of input reads, expressed as proportion (e.g., 0.03). If "0", then have flye set the read error automatically (default: 0)

Merge options:
  -I PERCENT, --circlator_min_id PERCENT
                        Percent identity threshold for circlator merge (default: 99)
  -L LENGTH, --circlator_min_length LENGTH
                        Minimum required overlap (bp) between original and merge contigs (default: 10000)
  -e LENGTH, --circlator_ref_end LENGTH
                        Minimum distance (bp) between end of original contig and nucmer hit (default: 100)
  -E LENGTH, --circlator_reassemble_end LENGTH
                        Minimum distance (bp) between end of merge contig and nucmer hit (default: 100)

Workflow options:
  -k, --keep_going_with_failed_contigs
                        Set this flag to continue running this script even if some contigs cannot be circularized and end-repaired. For any non-repaired contigs, an exact copy
                        of the original contig will be output.
  -c, --custom_assembly_info_file
                        Set this flag if you provide a custom tab-separated file to specify the circular vs. linear status of each contig in the assembly, rather than using
                        the assembly_info.txt file output by Flye (default), as your argument to -i. The custom file must have the following format: no headers; tab-separated;
                        first column is the contig names; second column is the status of the contigs, either "circular" or "linear". Any contig names not in this file will be
                        dropped by the script!
  -T LIST, --length_thresholds LIST
                        Comma-separated list of length thresholds for reassembly around the contig ends (bp) (default: 100000,75000,50000,25000,5000,2500,1000)
  -t JOBS, --threads JOBS
                        Number of processors threads to use (default: 1)
  -m GB, --threads_mem GB
                        Memory (GB) to use **per thread** for samtools sort (default: 1)
```

### `rotate`
Support module for rotating DNA sequences in FastA files.

Help menu (`spokewrench rotate -h`):
```commandline
usage: spokewrench rotate [-h] [-lf PATH] [-O] [-v] -i PATH -o PATH [-m] [-p INT] [-P FLOAT] [-t PATH] [-T PATH] [-n LIST] [-r PATH] [-s]

optional arguments:
  -h, --help            show this help message and exit

Basic config settings:
  -lf PATH, --logfile PATH
                        Log filepath (default: None)
  -O, --overwrite       Overwrite existing files/directories. By setting this flag, you risk erasing old data.
  -v, --verbose         Enable verbose logging

Required:
  -i PATH, --input_fasta PATH
                        Input fasta file
  -o PATH, --output_fasta PATH
                        Output fasta file. Note that all input sequences will always be written to output, even if only a subset of sequences are selected for rotate
                        operations.

Rotate options:
  -m, --midpoint        Rotate all sequences to their midpoint. (Incompatible with -p, -P, -t, and -T.)
  -p INT, --rotate_position INT
                        Number of base positions to rotate all sequences by, in the counter-clockwise direction. To specify a unique rotate position for each sequence, see
                        -t). Incompatible with -m, -P, -t, and -T.
  -P FLOAT, --rotate_fraction FLOAT
                        Fractional position to rotate all sequences to, counter-clockwise. For example, 0.3 means to rotate all contigs counter-clockwise to 30 percent of
                        their total length). To specify a unique fractional position for each sequence, see -T). Incompatible with -m, -p, -t, and -T.
  -t PATH, --rotate_position_table PATH
                        Path to a tab-separated file that contains the desired amount to rotate each sequence by, in base pair position numbers, in the counter-clockwise
                        direction. Columns (with headers) should be "sequence_id" and "rotate_bp". Only sequences specified in the table will be rotated, although all
                        sequences in the input file will be written to the output file. Incompatible with -m, -p, -P, and -T.
  -T PATH, --rotate_fraction_table PATH
                        Path to a tab-separated file that contains the desired amount to rotate each sequence by, as a fraction of total sequence length, in the counter-
                        clockwise direction. Columns (with headers) should be "sequence_id" and "rotate_fraction". Only sequences specified in the table will be rotated,
                        although all sequences in the input file will be written to the output file. Incompatible with -m, -p, -P, and -t.
  -n LIST, --sequence_names LIST
                        Sequences to be rotated (comma-separated list of IDs). Rotate operations will only be applied to these sequences, although all sequences in the input
                        file will be written to output. If used with -t or -T, the sequences listed in the table will be further subset to those that are also in the provided
                        list of sequence names.
  -r PATH, --output_report PATH
                        Path to write an optional tab-separated report file that shows how sequences were rotated.

Workflow options:
  -s, --strip_descriptions
                        Strip descriptions off FastA headers (i.e., any text after the first whitespace in each header)
```
