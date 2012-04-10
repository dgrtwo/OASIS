"""Constants.py
Created by David Robinson
6/26/09

Contains constants for OASIS program"""

### CONSTANTS ###

## BLAST executables
BLAST_EXE = "/Genomics/grid/users/dgrtwo/OASIS/bin/blastall"
FORMAT_EXE = "/Genomics/grid/users/dgrtwo/OASIS/bin/formatdb"

### ISFINDER LIBRARY FILES
ISFINDER_NUCLEOTIDE_FILE = "ISFinder_nuc_081311.fasta"
ISFINDER_AA_FILE = "ISFinder_aa_081311.fasta"

MINLENGTH = 600

LENGTH_WIGGLE = 2
MIN_DIVERGENCE = .95
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

## checking for inverted repeats around the edges of annotated IS elements
IR_WINDOW = 20
MIN_IR_SCORE_SMALL = 6

#identifying transposase family
TPASE_MAX_E_VALUE = 10 ** -12

MIN_IR_SCORE = 11

MIN_SCORE = 40
