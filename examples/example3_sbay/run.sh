#!/bin/bash
cd "$(dirname "$0")"
python3 ../../hybridchromomap.py -k karyotype.tsv -s segments.tsv -o output.png --width 14
