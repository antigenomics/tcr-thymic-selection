#!/bin/bash
set -e
set -x
SCRIPT_DIR="$(dirname ${BASH_SOURCE[0]})"
cd "${SCRIPT_DIR}/D1"
mkdir -p ../D1_intermediate_files
rm -f ../D1_intermediate_files/all_with_filenames.txt
for file in `ls -1 *.gz`; do
  zcat $file | sed -e "s/$/ $file/" >> ../D1_intermediate_files/all_with_filenames.txt
done
cd ../D1_intermediate_files
sort -u all_with_filenames.txt > all_without_duplicates.txt
cut -d" " -f1 all_without_duplicates.txt > all_filenames_stripped.txt
sort all_filenames_stripped.txt | uniq -c > counted_by_files.txt
cd ../D1
rm -f ../D1_intermediate_files/all.txt
zcat *.gz >> ../D1_intermediate_files/all.txt
cd ../D1_intermediate_files
sort all.txt | uniq -c > counted_by_lines.txt
join -j 2 -o 1.2,1.1,2.1 <(sort -k2 counted_by_files.txt) <(sort -k2 counted_by_lines.txt) > merged.txt
sort -k2,2 -k3,3 -r -n merged.txt | gzip -9 > ../sorted_d1.txt.gz
