#!/usr/bin/env bash
set -euo pipefail
# test-rotate-end-to-end.sh
# Runs end to end tests for the rotate module of spokewrench
# Copyright Jackson M. Tsuji, 2024

# If input field is not the appropriate length, print help and end script
if [ $# != 5 ]; then
  echo "$(basename $0): Run end to end test for the rotate module of spokewrench"
  echo ""
  echo "Usage: $(basename $0) input.fasta guide-positions.tsv guide-fractions.tsv expected_dir output_dir"
  echo ""
  echo "Positional arguments:"
  echo "   input.fasta: input FastA file for test sequence rotation"
  echo "   guide-positions.tsv: position-based guide file for test rotation"
  echo "   guide-fractions.tsv: fraction-based guide file for test rotation"
  echo "   expected_dir: directory containing expected output files from testing"
  echo "   output_dir: directory where output files will be saved"
  echo ""
  exit 0
fi

# Set input variables
input_fasta="$1"
guide_positions="$2"
guide_fractions="$3"
expected_dir="$4"
output_dir="$5"

# Test input variables
if [[ ! -f "${input_fasta}" ]]; then
  echo "ERROR: Input file '${input_fasta}' does not exist."
  exit 1
elif [[ ! -f "${guide_positions}" ]]; then
  echo "ERROR: Input file '${guide_positions}' does not exist."
  exit 1
elif [[ ! -f "${guide_fractions}" ]]; then
  echo "ERROR: Input file '${guide_fractions}' does not exist."
  exit 1
elif [[ ! -d "${expected_dir}" ]]; then
  echo "ERROR: Input directory '${expected_dir}' does not exist."
  exit 1
elif [[ -d "${output_dir}" ]]; then
  echo "ERROR: Output directory '${output_dir}' already exists."
  exit 1
fi

echo "Running tests..."
mkdir -p "${output_dir}"

# Do not exit early if errors are produced during testing
set +e

# Test series 1: midpoint rotation:
# Test 1a: midpoint rotate, keeping sequence comment lines (should finish)
spokewrench rotate -i "${input_fasta}" -o "${output_dir}/test1a.fasta" -r "${output_dir}/test1a.tsv" -v -m \
  2> "${output_dir}/test1a.log"; echo $? > "${output_dir}/test1a.exit"

# Test 1b: midpoint rotate, deleting sequence comment lines (should finish)
spokewrench rotate -i "${input_fasta}" -o "${output_dir}/test1b.fasta" -r "${output_dir}/test1b.tsv" -v -m -s \
  2> "${output_dir}/test1b.log"; echo $? > "${output_dir}/test1b.exit"

# Test 1c: target specific sequences (should finish)
spokewrench rotate -i "${input_fasta}" -o "${output_dir}/test1c.fasta" -r "${output_dir}/test1c.tsv" -v -m \
  -n "sequence2,sequence3" \
  2> "${output_dir}/test1c.log"; echo $? > "${output_dir}/test1c.exit"

# Test 1d: confirm that multiple rotation arguments makes the tool fail (should fail)
spokewrench rotate -i "${input_fasta}" -o "${output_dir}/test1d.fasta" -r "${output_dir}/test1d.tsv" -v -m -p 1 \
  2> "${output_dir}/test1d.log"; echo $? > "${output_dir}/test1d.exit"

# Test 1e: confirm the tool fails when an output FastA file is already present (should fail)
cp "${output_dir}/test1a.fasta" "${output_dir}/test1e.copied.fasta"
spokewrench rotate -i "${input_fasta}" -o "${output_dir}/test1e.copied.fasta" -r "${output_dir}/test1e.tsv" -v -m \
  2> "${output_dir}/test1e.log"; echo $? > "${output_dir}/test1e.exit"

# Test 1f: confirm the tool fails when an output TSV file is already present (should fail)
cp "${output_dir}/test1a.tsv" "${output_dir}/test1f.copied.tsv"
spokewrench rotate -i "${input_fasta}" -o "${output_dir}/test1f.fasta" -r "${output_dir}/test1f.copied.tsv" -v -m \
  2> "${output_dir}/test1f.log"; echo $? > "${output_dir}/test1f.exit"

# Test 1g: confirm the tool can continue when an output file is present but the overwrite flag is supplied (should finish)
echo "Example text" > "${output_dir}/test1g.fasta"
spokewrench rotate -i "${input_fasta}" -o "${output_dir}/test1g.fasta" -r "${output_dir}/test1g.tsv" -v -m -O \
  2> "${output_dir}/test1g.log"; echo $? > "${output_dir}/test1g.exit"


# Test series 2: rotate by common position
# Test 2a: normal rotation by position (should finish)
spokewrench rotate -i "${input_fasta}" -o "${output_dir}/test2a.fasta" -r "${output_dir}/test2a.tsv" -v -p 1 \
  2> "${output_dir}/test2a.log"; echo $? > "${output_dir}/test2a.exit"

# Test 2b: rotate specific sequences (should finish)
spokewrench rotate -i "${input_fasta}" -o "${output_dir}/test2b.fasta" -r "${output_dir}/test2b.tsv" -v -p 1 \
  -n "sequence2,sequence3" \
  2> "${output_dir}/test2b.log"; echo $? > "${output_dir}/test2b.exit"

# Test 2c: rotate to position 0 (should not rotate) (should finish)
spokewrench rotate -i "${input_fasta}" -o "${output_dir}/test2c.fasta" -r "${output_dir}/test2c.tsv" -v -p 0 \
  2> "${output_dir}/test2c.log"; echo $? > "${output_dir}/test2c.exit"

# Test 2d: rotate to negative position (should fail)
spokewrench rotate -i "${input_fasta}" -o "${output_dir}/test2d.fasta" -r "${output_dir}/test2d.tsv" -v -p -1 \
  2> "${output_dir}/test2d.log"; echo $? > "${output_dir}/test2d.exit"

# Test 2e: rotate to position larger than all sequences (should fail)
spokewrench rotate -i "${input_fasta}" -o "${output_dir}/test2e.fasta" -r "${output_dir}/test2e.tsv" -v -p 100 \
  2> "${output_dir}/test2e.log"; echo $? > "${output_dir}/test2e.exit"

# Test 2f: rotate to position larger than some sequences (should fail)
spokewrench rotate -i "${input_fasta}" -o "${output_dir}/test2f.fasta" -r "${output_dir}/test2f.tsv" -v -p 10 \
  2> "${output_dir}/test2f.log"; echo $? > "${output_dir}/test2f.exit"

# Test 2g: rotate to the same length as the shortest contig (should finish)
spokewrench rotate -i "${input_fasta}" -o "${output_dir}/test2g.fasta" -r "${output_dir}/test2g.tsv" -v -p 8 \
  2> "${output_dir}/test2g.log"; echo $? > "${output_dir}/test2g.exit"


# Test series 3: rotate by common fraction
# Test 3a: normal rotation by fraction (should finish)
spokewrench rotate -i "${input_fasta}" -o "${output_dir}/test3a.fasta" -r "${output_dir}/test3a.tsv" -v -P 0.3 \
  2> "${output_dir}/test3a.log"; echo $? > "${output_dir}/test3a.exit"
# Most other cases are covered in the Test 2 series.


# Test series 4: rotate by sequence-specific info via a guide file
# Test 4a: rotate by position via a guide file (should finish)
spokewrench rotate -i "${input_fasta}" -o "${output_dir}/test4a.fasta" -r "${output_dir}/test4a.tsv" -v \
  -t "${guide_positions}" \
  2> "${output_dir}/test4a.log"; echo $? > "${output_dir}/test4a.exit"
# TODO - consider changing the column names of pos.tsv and confirmed an error is produced

# Test 4b: rotate by position via a guide file but only for selected sequences (should finish)
spokewrench rotate -i "${input_fasta}" -o "${output_dir}/test4b.fasta" -r "${output_dir}/test4b.tsv" -v \
  -t "${guide_positions}" -n "sequence2,sequence3" \
  2> "${output_dir}/test4b.log"; echo $? > "${output_dir}/test4b.exit"

# Test 4c: rotate by fraction via a guide file (should finish)
spokewrench rotate -i "${input_fasta}" -o "${output_dir}/test4c.fasta" -r "${output_dir}/test4c.tsv" -v \
  -T "${guide_fractions}" \
  2> "${output_dir}/test4c.log"; echo $? > "${output_dir}/test4c.exit"

set -e

# See if all expected output files were generated
find "${output_dir}" -name "*.fasta" -or -name "*.tsv" -or -name "*.exit" | grep -v "copied" | \
  sed -e "s|^${output_dir}/||g" | sort -h > "${output_dir}/outputs.list"
find "${expected_dir}" -name "*.fasta" -or -name "*.tsv" -or -name "*.exit" | grep -v "copied" | \
  sed -e "s|^${expected_dir}/||g" | sort -h > "${output_dir}/expected.list"
set +e
cmp_status=$(cmp "${output_dir}/outputs.list" "${output_dir}/expected.list" >/dev/null 2>&1 && echo $?)
set -e
if [[ "${cmp_status}" != 0 ]]; then
  # TODO: If any single run fails, some files like the output .fasta will not be generated, causing an error.
  #       Need to make the testing more specific.
  echo "ERROR: output files from the test do not match the expected output files in '${expected_dir}'."
  echo "Please compare '${output_dir}/outputs.list' to '${output_dir}/expected.list' for troubleshooting."
  exit 1
fi

# Check the contents of each output file
failed_tests=0
test_files=($(cat "${output_dir}/outputs.list"))
for test_filename in "${test_files[@]}"; do
  test_filepath="${output_dir}/${test_filename}"
  ref_filepath="${expected_dir}/${test_filename}"

  set +e
  cmp_status=$(cmp "${ref_filepath}" "${test_filepath}" >/dev/null 2>&1 && echo $?)
  set -e

  if [[ "${cmp_status}" == 0 ]]; then
    echo "[ PASS ]: ${test_filename}"
  else
    echo "[ FAIL ]: ${test_filename}"
    failed_tests=$((failed_tests+1))
  fi
done

# Status report
echo "###########################"
if [[ $failed_tests == 0 ]]; then
  echo "SUMMARY: All tests passed."
  rm -r "${output_dir}"
  echo "###########################"
  echo "$(basename $0): Done."
else
  echo "SUMMARY: ${failed_tests} tests failed. Will not remove output dir '${output_dir}'."
  echo "###########################"
  echo "$(basename $0): Done, but exiting with status code 1 due to test failures."
  exit 1
fi
