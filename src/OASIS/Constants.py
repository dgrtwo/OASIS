"""Constants.py
Created by David Robinson
6/26/09

Contains constants for OASIS program"""

### CONSTANTS ###

## BLAST executables
import os
import ConfigParser

config_file = os.path.join(os.path.split(__file__)[0], "data", "data.cfg")
config_parser = ConfigParser.ConfigParser()
config_parser.read(config_file)

BLAST_EXE = config_parser.get("BLAST", "BLAST_EXE")
FORMAT_EXE = config_parser.get("BLAST", "FORMAT_EXE")

if BLAST_EXE.startswith("/PATH/TO") or FORMAT_EXE.startswith("/PATH/TO"):
    raise Exception("Must set full path to NCBI blastall executable in " +
                    config_file)

ISFINDER_AA_FILE = "ISFinder_aa.fasta"

# temporary directory
TEMPORARY_DIRECTORY = "temp_OASIS"

if not os.path.exists(TEMPORARY_DIRECTORY):
    os.mkdir(TEMPORARY_DIRECTORY)


MINLENGTH = 600

LENGTH_WIGGLE = 2
MIN_AA_SIMILARITY = .95
EXTEND_LENGTH = 25
MAX_MISS = 3
FULL_EXTEND = 2000
FULL_MAX_MISS = 2
REDUNDANT_WIGGLE = 20

FEATURE_WIGGLE = 10

## maximum window to consider two adjacent coding regions ORFs of the same gene
ORF_WINDOW = 100

## Maximum overlap to count two IS elements as the same, for re-annotating groups
MAX_OVERLAP = 300

## Constants used in BLAST for partial sequences
E_VALUE_CUTOFF = pow(10, -30) # maximum e value for a match to count as a partial
MIN_PARTIAL_LEN = 250 # minimum length for a match to count as a partial

## windows to check for inverted repeats for single copies:
SINGLE_IR_OUT_WINDOW = 500
SINGLE_IR_IN_WINDOW = 20
SINGLE_IR_MIN_SCORE = 11

## checking for inverted repeats around the edges of annotated IS elements
IR_WINDOW = 100
MIN_IR_SCORE_CHANGE = 11

# checking for inverted repeats immediately on the inside- lower penalty
IR_WINDOW_NONCHANGE = 20
MIN_IR_SCORE_NONCHANGE = 10


## families that do not have inverted repeats
NON_IR_FAMILIES = ["IS110", "IS200/IS605"]

#identifying transposase family
TPASE_MAX_E_VALUE = 10 ** -12

MIN_IR_SCORE = 11

MIN_SCORE = 40
