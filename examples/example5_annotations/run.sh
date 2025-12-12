#!/bin/bash
cd "$(dirname "$0")"
# With labels
python3 ../../hybridchromomap.py -k karyotype.tsv -s segments.tsv -c colors.tsv -a annotations.tsv -o output_with_labels.png
# Without labels
python3 ../../hybridchromomap.py -k karyotype.tsv -s segments.tsv -c colors.tsv -a annotations.tsv -o output.png --no-labels
