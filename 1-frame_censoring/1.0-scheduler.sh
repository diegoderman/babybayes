#!/bin/bash

###############################
#
# system-dependent scheduler configuration goes here
#
###############################

python3 ./calculate-dvars.py --in ../data/table.csv --out ../data/table_dvars.csv --threads 48
python3 ./filter-dvars.py --in ../data/table.csv --out ../data/table_dvars.csv --threads 48
