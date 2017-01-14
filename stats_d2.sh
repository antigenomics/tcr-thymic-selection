#!/bin/bash
set -e
set -x
SCRIPT_DIR="$(dirname ${BASH_SOURCE[0]})"
cd "${SCRIPT_DIR}/D2"
mkdir -p ../D2_intermediate_files
rm -f ../D2_intermediate_files/cdr3aa.txt
rm -f ../D2_intermediate_files/cdr3aa_with_filenames.txt
for file in `ls -1 *.gz`; do
  zcat $file | tail -n +2 | cut -f4 >> ../D2_intermediate_files/cdr3aa.txt
  zcat $file | tail -n +2 | cut -f4 | sed -e "s/$/ $file/" >> ../D2_intermediate_files/cdr3aa_with_filenames.txt
done
cd ../D2_intermediate_files
sort -u cdr3aa_with_filenames.txt > cdr3aa_without_duplicates.txt
cut -d" " -f1 cdr3aa_without_duplicates.txt > cdr3aa_filenames_stripped.txt
sort cdr3aa_filenames_stripped.txt | uniq -c > counted_by_files.txt
cd ../D2_intermediate_files
sort cdr3aa.txt | uniq -c > counted_by_lines.txt
join -j 2 -o 1.2,1.1,2.1 <(sort -k2 counted_by_files.txt) <(sort -k2 counted_by_lines.txt) > merged.txt
sort -k2,2 -k3,3 -r -n merged.txt | gzip -9 > ../sorted_d2.txt.gz
