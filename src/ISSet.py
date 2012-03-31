"""ISSet.py
Created by David Robinson
8/11/08

OASIS is a module designed for IS element annotation for prokaryote 
genomes.

This module represents a set of copies of an IS element."""

#imports
import os

from Bio import SeqRecord

import IS
import EdgeFinder
from OASIS_functions import *
from Constants import *

class ISSet:
    """represents a set of IS matches that are one type"""
    def __init__(self, lst, profile):
        """initialized with a list of IS objects and a profile object"""
        self.lst = lst
        self.profile = profile
        self.type = None
        first = self.lst[0]
        if first.feature and 'product' in first.feature.qualifiers:
            product = first.feature.qualifiers['product'][0]
            if "transposase of " in product:
                self.type = product.split("sase of ")[1].split(",")[0]
    
    def re_annotate(self, from_edges=False, window=FULL_EXTEND):
        """create a position specific scoring matrix for the surrounding regions,
        see how far you can go while keeping conservation, using a maximum
        likelihood method. Can either annotate from feature (by default) or
        from edges"""
        if len(self.lst) < 2:
            return
        
        #unzip
        if from_edges:
            [befores, afters] = zip(*[e.around_IS(window) for e in self.lst])
        else:
            [befores, afters] = zip(*[e.around_gene(window) for e in self.lst])
        
        # make them the same length
        minlength_before = min([len(b) for b in befores])
        minlength_after = min([len(a) for a in afters])
        
        befores = [b[:minlength_before] for b in befores]
        afters = [a[:minlength_after] for a in afters]
        
        # old:
        #before_len = multiple_extend(befores)
        #after_len = multiple_extend(afters)
        before_len = EdgeFinder.find_left_edge([str(s) for s in befores])
        after_len = EdgeFinder.find_left_edge([str(s) for s in afters])
        
        if from_edges and False:
            print "BEFORES:"
            print "\n".join([str(b[:20]) + "\t" + str(e) for b, e in zip(befores, self.lst)])
            print "AFTERS:"
            print "\n".join([str(a[:20]) + "\t" + str(e) for a, e in zip(afters, self.lst)])
            print "before_len", before_len
            print "after_len", after_len
        
        # TESTING
        if False:
            print before_len
            print after_len
            
            for s in befores:
                print str(s[before_len-20:before_len]) + "|" + str(s[before_len:before_len+20])
            
            for s in afters:
                print str(s[after_len-20:after_len]) + "|" + str(s[after_len:after_len+20])
            
            print
            print
        
        for e in self.lst:
            if from_edges:
                start = e.start
                end = e.end+1
            else:
                start = e.feature.location._start.position
                end = e.feature.location._end.position
            
            if e.direction == 1:
                newstart = start - before_len
                newend = end + after_len
            else:
                newstart = start - after_len
                newend = end + before_len
            
            e.change(newstart=newstart, newend=newend)
                
    def add_IRs(self):
        """add inverted repeats to each IS in this set"""
        [e.set_IRs() for e in self.not_partial()]
    
    def filter(self):
        """include IS elements that meet the conditions- return False if the
        overall group does not meet the conditions"""
        self.lst = [e for e in self.lst if e.meets_conditions()]
        if len(self.lst) >= 4: return True #is very rarely a mistake
        return len(self.lst) > 0 and self.lst[0].is_valid()
    
    def add(self, other_type):
        """combines the other set of IS elements with this one"""
        self.lst = self.lst + other_type.lst
    
    def not_partial(self):
        """return a list of all IS elements that are not partials"""
        return [e for e in self.lst if e.length >= MINLENGTH]
    
    def representative(self):
        """return a representative IS element for blasting"""
        choices = self.lst
        for e in choices:
            if e.IRL and len(e.IRL) > 8:
                choices = [e for e in choices if e.IRL and len(e.IRL) > 8]
                break
        lens = [e.length for e in choices]
        lsorted = list(set(lens))
        lsorted.sort(lambda a, b: cmp(lens.count(b), lens.count(a)))
        return choices[lens.index(lsorted[0])]
    
    def remove_redundant(self, genome):
        """remove IS elements that are already annotated in the genome. Returns
        whether there are any IS elements left"""
        self.lst = [e for e in self.lst if not genome.redundant_IS(e)]
        return self.lst > 0
    
    def score(self):
        """return a score representing how good a group this is"""
        score = 0
        for e in self.lst:
            if e.IRL: score = score + len(e.IRL)
        return score
    
    def identify_transposase(self):
        """set the families and groups of each or False if there is no 
        family or group"""
        aaseqs = [e.aaseq for e in self.lst if e.aaseq]
                
        aaseqs.sort(lambda a, b: cmp(len(b), len(a)))
        if len(aaseqs) == 0: return False
        
        family, group = self.profile.identify_family(aaseqs[0])
        
        #set variables
        for e in self.lst: e.family = family
        for e in self.lst: e.group = group
        
        return
        
        #create temporary fasta file
        tpase_db = "ISFinder_transposases.fasta"
        blast_file = "temp.fasta"
        outf = open(blast_file, "w")
        
        id = "temp"
        aaseqs = [e.aaseq for e in self.lst if e.aaseq]
        aaseqs.sort(lambda a, b: cmp(len(b), len(a)))
        if len(aaseqs) == 0: return False
        
        temp_record = SeqRecord.SeqRecord(id=id, seq=aaseqs[0])
        
        SeqIO.write([temp_record], outf, "fasta")
        outf.close()
        
        #perform blast
        result_handle, error_handle = NCBIStandalone.blastall(BLAST_EXE,
                                        "blastp", tpase_db, blast_file)
        
        record = NCBIXML.parse(result_handle).next()
        
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
            fields = best_alignment.title.split(" ")[1].split("|")
            if fields[2] != "" and fields[2] != "-":
                family = fields[2]
            if fields[3] != "" and fields[3] != "-":
                group = fields[3]
                
        #clean up by removing temporary blast file
        os.system("rm " + blast_file)
    
    def fasta(self):
        """return as a list of SeqRecord objects, one for each IS"""
        retlist = []
        for e in self.lst:
            retlist = retlist + [e.fasta()] + e.aa_fasta()
        return retlist
        
    def __repr__(self):
        """return all match strings, ended with a newline"""
        return "".join ([str(m) for m in self.lst]) + "\n"