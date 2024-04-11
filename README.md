# rotary-utils
Circularization-related utilities for the rotary pipeline

## Overview
*rotary-utils*, currently under development, is a suite of utilities for manipulating circular DNA elements.
Once completed, these utilities will be able to serve as a modern replacement for
[circlator](https://github.com/sanger-pathogens/circlator), which is now in a frozen development state.
The utilities in *rotary-utils* are built into the [*rotary* project](https://github.com/rotary-genomics/rotary),
where they are used with "best practices" to ensure proper assembly of circular sequences. These utilities can also
be used in a standalone fashion or integrated into custom genome assembly workflows to ensure circular DNA/RNA elements
are assembled accurately, similar to how circlator is currently used.

## Requirements
- OS: Runs on Linux (tested on Ubuntu 20.04 and Ubuntu 22.04) and macOS (tested on Sonoma 14)
- Software: requires miniconda or manual installation using the dependencies shown in rotary-utils/environment.yml
- Resources: should run on a modern laptop with >=8 GB RAM and >=4 CPU threads, in most cases

## Installation
```commandline
git clone https://github.com/rotary-genomics/rotary-utils.git

conda env create -n rotary_utils --file=rotary-utils/environment.yml

conda activate rotary_utils

cd rotary-utils

pip install --editable .

# See available commands
rotary-utils -h
```

## Modules
### `repair`
Run on the outputs of Flye to repair the
[short gap or overlap region](https://github.com/fenderglass/Flye/issues/315#issuecomment-720679812) that can occur at
the ends of circular contigs produced by this assembly tool.

Help menu (`rotary-utils repair -h`):
```commandline
usage: rotary-utils repair [-h] -l LONG_READ_FILEPATH -a ASSEMBLY_FASTA_FILEPATH -i ASSEMBLY_INFO_FILEPATH -o OUTPUT_DIR [-f {nano-hq,nano-raw}] [-F FLYE_READ_ERROR]
                           [-I CIRCLATOR_MIN_ID] [-L CIRCLATOR_MIN_LENGTH] [-e CIRCLATOR_REF_END] [-E CIRCLATOR_REASSEMBLE_END] [-k] [-c] [-O] [-T LENGTH_THRESHOLDS]
                           [-t THREADS] [-m THREADS_MEM] [-v]

optional arguments:
  -h, --help            show this help message and exit

Required:
  -l LONG_READ_FILEPATH, --long_read_filepath LONG_READ_FILEPATH
                        QC-passing Nanopore reads
  -a ASSEMBLY_FASTA_FILEPATH, --assembly_fasta_filepath ASSEMBLY_FASTA_FILEPATH
                        Contigs to be end-repaired
  -i ASSEMBLY_INFO_FILEPATH, --assembly_info_filepath ASSEMBLY_INFO_FILEPATH
                        assembly_info.txt file from Flye showing which assembled contigs are circular vs. linear (or a custom guide file; see -c below)
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory path

Input read options:
  -f {nano-hq,nano-raw}, --flye_read_mode {nano-hq,nano-raw}
                        Type of long read reads provided by -l, to be used for reassembly by Flye. See details on these settings in the Flye documentation. (default: nano-hq)
  -F FLYE_READ_ERROR, --flye_read_error FLYE_READ_ERROR
                        Expected error rate of input reads, expressed as proportion (e.g., 0.03). If "0", then have flye set the read error automatically (default: 0)

Merge options:
  -I CIRCLATOR_MIN_ID, --circlator_min_id CIRCLATOR_MIN_ID
                        Percent identity threshold for circlator merge (default: 99)
  -L CIRCLATOR_MIN_LENGTH, --circlator_min_length CIRCLATOR_MIN_LENGTH
                        Minimum required overlap (bp) between original and merge contigs (default: 10000)
  -e CIRCLATOR_REF_END, --circlator_ref_end CIRCLATOR_REF_END
                        Minimum distance (bp) between end of original contig and nucmer hit (default: 100)
  -E CIRCLATOR_REASSEMBLE_END, --circlator_reassemble_end CIRCLATOR_REASSEMBLE_END
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
  -O, --overwrite       Set this flag to overwrite an existing output directory and its contents; by default this code will throw an error if the output directory already
                        exists. By setting this flag, you risk erasing old data and/or running into errors.
  -T LENGTH_THRESHOLDS, --length_thresholds LENGTH_THRESHOLDS
                        Comma-separated list of length thresholds for reassembly around the contig ends (bp) (default: 100000,75000,50000,25000,5000,2500,1000)
  -t THREADS, --threads THREADS
                        Number of processors threads to use (default: 1)
  -m THREADS_MEM, --threads_mem THREADS_MEM
                        Memory (GB) to use **per thread** for samtools sort (default: 1)
  -v, --verbose         Enable verbose logging
```
