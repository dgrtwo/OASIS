"""OASIS_functions.py
Created by David Robinson
8/11/08

OASIS is a module designed for IS element annotation for prokaryote
genomes.

This module creates useful functions that are used by various modules."""

import os
import re
import math
import random

from Bio import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC

from OASIS_functions import *
from Constants import *
from SW import myalign
import my_SW

#set up parameters
nucs = ["A", "G", "C", "T"]

#functions
def related_text(txt, single):
    """return whether a protein is a transposase or integrase from its
    product description"""
    return (("transposase" in txt or ("integrase" in txt and not single))
                and ("integrase family" not in txt))

def related(feature, single=False):
    """recognizes whether a feature is relevant"""
    if feature.type != "CDS":
        return False
    if "product" in feature.qualifiers:
        return related_text(feature.qualifiers["product"][0].lower(), single)
    if "function" in feature.qualifiers:
        return related_text(feature.qualifiers["function"][0].lower(), single)
    return False

def extend(seq1, seq2, maxmis):
    """given two sequences and the maximum mismatches permitted for an
    alignment, return the length of alignment that can be achieved"""
    i = 0
    mm = 0
    lastgood = 0 #keep track of the last valid bp
    while i < len(seq1) and i < len(seq2) and mm <= maxmis:
        if seq1[i] != seq2[i]:
            mm = mm + 1
        elif seq1[i-1] == seq2[i-1]:
            lastgood = i
        i = i + 1
    if lastgood == len(seq1)-1:
        return lastgood + 1
    return lastgood

def IS_in_list(m, lst):
    """check whether this match is close to any matches in the list"""
    for e in lst:
        if m.chromosome == e.chromosome and \
            m.start >= e.start - REDUNDANT_WIGGLE \
            and m.end <= e.end + REDUNDANT_WIGGLE:
            return True
    return False

def IRs(seq, verbose=False):
    """takes a Seq object representing an IS element sequence.
    Returns IR elements on either side, as well as their starting indices,
    and a similarity score"""
    start = str(seq[:IR_WINDOW])
    end = str(seq[-IR_WINDOW:].reverse_complement())

    #aln = pairwise2.align.localms(start, end, 1, -20, -5, -2)
    aln = myalign(start, end)

    if (aln[2] < MIN_IR_SCORE_CHANGE):
        # try a close alignment with a lower penalty- one that doesn't move
        # based on the alignment, and accepts only an exact match
        close_aln = myalign(start[:IR_WINDOW_NONCHANGE],
                            end[:IR_WINDOW_NONCHANGE], mismatch_score_num=-1)

        if (close_aln[2] < MIN_IR_SCORE_NONCHANGE or
            close_index(start, close_aln[0]) != 0 or
            close_index(end, close_aln[1]) != 0):
            # no alignment near or far
            return False, False, 0, 0, 0
        return close_aln[0], close_aln[1], 0, 0, close_aln[2]

    lin, rin = close_index(start, aln[0]), -close_index(end, aln[1])

    return aln[0], aln[1], lin, rin, aln[2]


def seq_match(seq1, seq2, threshold, mm):
    """return whether seq1 and seq2 align to threshold length)"""
    return extend(seq1[:threshold], seq2[:threshold], mm) == threshold


def close_index(string, substring):
    if substring in string: return string.index(substring)
    #otherwise, try successively shorter strings
    i = len(substring)
    while i > 0:
        this_substring = substring[:i]
        if this_substring in string: return string.index(this_substring)
        i = i - 1
    #failed- none of the common string is in this
    return -1


def similar(set1, set2, verbose=False):
    """given two is_sets, return whether they are probably the same IS element"""
    full_list1 = [e for e in set1.lst if e.length >= MIN_PARTIAL_LEN]
    full_list2 = [e for e in set2.lst if e.length >= MIN_PARTIAL_LEN]

    lsts = full_list1 + full_list2
    chromosomes = list(set([e.chromosome for e in lsts]))

    overlaps = 0
    all_indices = []

    for c in chromosomes:
        this_lst = [e for e in lsts if e.chromosome == c]
        this_lst.sort(key=lambda x: x.start)

        #print this_lst

        indices = [this_lst[i].end-this_lst[i+1].start for i in range(len(this_lst)-1)]

        all_indices = all_indices + indices

        if verbose:
            print [e for e in indices if e > MAX_OVERLAP]

        overlaps = overlaps + len([e for e in indices if e > MAX_OVERLAP])

    #print all_indices, overlaps

    #print [lsts[i+1].end-lsts[i].start for i in range(len(lsts)-1)]

#     if len(full_list1) != len(full_list2): return False
#     similar_count = 0
#     for e in full_list1:
#         if e.genename:
#             if e.genename in [e2.genename for e2 in full_list2]:
#                 similar_count = similar_count + 1
#         else:
#             if e.start in [e2.start for e2 in full_list2]:
#                 similar_count = similar_count + 1

    #is_similar = (float(similar_count) / float(len(full_list1)) > .1) or \

    if verbose:
        print overlaps

    fraction_similar = 10

    return (overlaps > len(full_list1)/fraction_similar or overlaps > len(full_list2)/fraction_similar)


def translate(seq):
    """return a sequence translated with the bacterial code"""
    newseq = Seq.Seq(str(seq), alphabet=IUPAC.unambiguous_dna)
    return newseq.translate()
    #return bacteria_translator.translate(newseq)


def same_gene(f1, f2):
    """check if two ORFs are in the same gene"""
    end1 = min([f.location._end.position for f in [f1, f2]])
    start2 = max([f.location._start.position for f in [f1, f2]])

    return start2 - end1 <= ORF_WINDOW


def similarity(seq1, seq2):
    """
    takes two sequences as strings and calculates their similarity with
    a simple algorithm
    """
    matchnum = 0
    i = 0
    j = 0
    while True:
        if seq1[i] == seq2[j]: matchnum = matchnum + 1
        else:
            #check for skip:
            for change in [3]:
                if seq1[i:i+change] == seq2[j+change:j+change+change]:
                    j = j + change - 1
                    i = i - 1
                if seq2[j:j+change] == seq1[i+change:i+change+change]:
                    i = i + change - 1
                    j = j - 1
        i = i + 1
        j = j + 1

        if i >= len(seq1) or j >= len(seq2): break
        if i >= 6 and matchnum < i/2: break

    return float(matchnum) / float(len(seq1))


def find_IRs(family, seq1, seq2, in_window):
    """given the sequence and family of an IS and the windows around it,
    find the most likely inverted repeats"""
    #change to strings
    window1 = str(seq1)
    window2 = str(seq2)

    start_i, max_i, start_j, max_j, score = my_SW.align(window1, window2)

    IR1 = seq1[start_i-1:max_i]
    IR2 = seq2[start_j-1:max_j].reverse_complement()

    if score > SINGLE_IR_MIN_SCORE:
        # return actual IR sequences
        return IR1, IR2
        return max_i-in_window, len(seq2)-max_j-in_window
    else:
        return False
