"""Profile.py
Created by David Robinson
10/22/08

Represents a profile of ISFinder IS families"""

#imports
import os
import re

from Bio import SeqRecord
from Bio.Blast import NCBIStandalone
from Bio.Blast import NCBIXML
from Bio import SeqIO

from OASIS_functions import *
from Constants import *
import my_SW

#classes

class Profile:
    """
    describes a profile of ISFinder transposases that can identify their
    gene
    """
    def __init__(self, tpase_file):
        """initialized with a sequence file and a transposase file"""
        #save file for future reference
        self.tpase_file = tpase_file
        
        #turn each into a database
        os.system(FORMAT_EXE + " -p T -i " + tpase_file)
    
    def identify_family(self, aaseq):
        """given an amino acid sequence, identify its family"""
        blast_file = "profile_temp.fasta"
        outf = open(blast_file, "w")
                
        temp_record = SeqRecord.SeqRecord(id="temp", seq=aaseq)
        
        SeqIO.write([temp_record], outf, "fasta")
        outf.close()
        
        #perform blast
        print BLAST_EXE, self.tpase_file, blast_file
        result_handle, error_handle = NCBIStandalone.blastall(BLAST_EXE,
                                        "blastp", self.tpase_file, blast_file)

        try:
            record = NCBIXML.parse(result_handle).next()
        except ValueError:
            print "NCBI BLAST error: " + error_handle.read()
            raise Exception()
        
        best_hsp = None
        best_alignment = None
        
        #perform blast
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < TPASE_MAX_E_VALUE:
                    if best_hsp:
                        if hsp.score > best_hsp.score:
                            best_alignment = alignment
                            best_hsp = hsp
                    else:
                        best_alignment = alignment
                        best_hsp = hsp
        
        #find family and group
        family = None
        group = None
        
        if best_hsp:
            fields = re.split("[\s\t]+", best_alignment.title)[1].split("|")
            #best_IS = self.__fetch_by_name(fields[0])
            family, group = fields[2], fields[3]
        
        #clean up by removing temporary blast file
        os.system("rm " + blast_file)
                
        return family, group
    
    def find_IRs(self, family, seq1, seq2, in_window):
        """given the sequence and family of an IS and the windows around it,
        find the most likely inverted repeats"""    
        #change to strings
        window1 = str(seq1)
        window2 = str(seq2)
        
        #print window1, window2
        
        start_i, max_i, start_j, max_j, score = my_SW.align(window1, window2)
        #print start_i, max_i, start_j, max_j, score
        
        #print score
        #print seq1[start_i-1:max_i]
        #print seq2[start_j-1:max_j]
                
        IR1 = seq1[start_i-1:max_i]
        IR2 = seq2[start_j-1:max_j].reverse_complement()
        
        if score > MIN_IR_SCORE_SMALL:
            # return actual IR sequences
            return IR1, IR2
            return max_i-in_window, len(seq2)-max_j-in_window
        else:
            return False


