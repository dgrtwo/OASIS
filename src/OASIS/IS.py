"""IS.py
Created by David Robinson
8/11/08

OASIS is a module designed for IS element annotation for prokaryote 
genomes.

This module represents the data on a genome that allows it to be annotated
This data can be imported from a genbank file or given to it as a list of
SeqRecords."""

#imports
from Bio import SeqIO
from Bio import SeqRecord

from OASIS_functions import *
from Constants import *

### PARAMETERS ###

DIRECTION_TO_STRAND = {1: "+", -1: "-"}

class IS:
    """represents an annotated IS element and some information about it"""
    def __init__(self, feature, chromosome, start, end, gen, 
                family=None, group=None, dir=None):
        self.feature = feature
        if not feature:
            #gene at this position is not known
            if not dir:
                raise Exception("IS must be provided with direction or a gene")
            self.genename = ""
            self.direction = dir
            self.aaseq = None
        else:
            self.genename = self.feature.qualifiers['locus_tag'][0]
            self.direction = {None: 1, 1: 1, -1:-1}[feature.strand]
            self.aaseq = gen.get_aaseq(self.feature)
        self.chromosome = chromosome
        self.start = start
        self.end = end
        if gen:
            self.chrom_seq = gen.get_chrom(chromosome)
            self.__setseq()
        
        self.family = family
        self.group = group
        
        if gen:
            allgenes = gen.get_multiple_features(chromosome, start, end)
            self.names = [f.qualifiers['locus_tag'][0] for f in allgenes]
            
            self.aaseqs = [gen.get_aaseq(f) for f in 
                            gen.get_multiple_features(chromosome, start, end)]
        
        #if self.genename == "Aave_0513":
        #    print self.genename, chromosome, start, end
        #    print map(str, gen.get_multiple_features(chromosome, start, end))
        
        #now some optional variables that can be set later:
        self.IRL = None
        self.IRR = None
    
    def __setseq(self):
        """set the self.seq and self.length attributes"""
        self.length = self.end - self.start + 1
        self.seq = self.chrom_seq[self.start:self.end]
    
    def change(self, startchange=0, endchange=0, absolute=False, \
               newstart=None, newend=None):
        """change an IS element's location and its sequence. If absolute, use
        the start and end based on IS direction, rather than which is smaller"""
        if newstart: self.start = newstart
        if newend: self.end = newend
        if not absolute or self.direction == 1:
            self.start = self.start + startchange
            self.end = self.end + endchange
        else:
            self.start = self.start - startchange
            self.end = self.end - endchange
        self.__setseq()
    
    def set_IRs(self):
        """find the inverted repeats, if any, around an IS"""
        if self.family in NON_IR_FAMILIES:
            return

        verbose = False
        s = self.seq if self.direction == 1 else self.seq.reverse_complement()
        self.IRL, self.IRR, lin, rin, self.IRscore = IRs(s, verbose)
        if verbose:
            print s
            print self.as_gff(1), self.IRL, self.IRR, lin, rin

        # if this is reversed, switch the adjustments around
        if self.direction == -1:
            lin, rin = -rin, -lin

        self.change(startchange=lin, endchange=rin, absolute=False)
        
    def meets_conditions(self):
        """return whether this IS is satisfactory"""
        return self.length >= MINLENGTH
    
    def is_valid(self):
        """check if this element appears to be a valid IS element"""
        #first, check if it has the word "transposase"
        if "product" in self.feature.qualifiers:
            product = self.feature.qualifiers['product'][0]
            if "transposase" in product or "integrase" in product:
                return True
        #otherwise, see if it has valid IRs
        if self.IRL and len(self.IRL) > 6: return True
        #finally, try blasting
        #blast_tpase = is_transposase(self.aaseq)
        #if not blast_tpase:
        #    print self.feature.qualifiers['product'][0]
        return False
    
    def around_gene(self, window):
        gene_start = self.feature.location._start.position
        gene_end = self.feature.location._end.position
        before = self.chrom_seq[max([gene_start-window,0]):gene_start].reverse_complement()
        after = self.chrom_seq[gene_end:gene_end+window]
        if self.direction == -1:
            before, after = after, before
        return before, after
    
    def around_IS(self, window):
        start = self.start
        end = self.end
        before = self.chrom_seq[max([start-window,0]):start].reverse_complement()
        after = self.chrom_seq[end:end+window]
        if self.direction == -1:
            before, after = after, before
        return before, after

    def __repr__(self):
        retstr = ""
        #add basic information on chromosome, start, end, length, genome name
        retstr = retstr + self.chromosome + "\t"
        retstr = retstr + str(self.start+1) + "\t" + str(self.end) + "\t"
        retstr = retstr + str(self.end-self.start) + "\t" + str(self.direction)
        
        retstr = retstr + "\t" + "/".join(self.names) + "\t"
        
        if self.feature and'product' in self.feature.qualifiers:
            product = self.feature.qualifiers['product'][0]
            istype = ""
            if "transposase of " in product:
                istype = product.split("sase of ")[1].split(",")[0]
            retstr = retstr + product + "\t" + istype + "\t"
        else:
            retstr = retstr + "\t\t"
        
        retstr = retstr + str(self.family) + "\t" + str(self.group) + "\t"
        
        if self.IRL:
            retstr = retstr + self.IRL + "\t" + self.IRR + "\t"
        else:
            retstr = retstr + "\t\t"
        
        # temporarily add surrounding regions
        #before = self.chrom_seq[self.start-20:self.start]
        #after = self.chrom_seq[self.end:self.end+20]
        #retstr = retstr + str(before) + "\t" + str(after) + "\t"
        
        retstr = retstr + "\n"
        #retstr = retstr + str(self.upstream_len) + "\t" + str(self.downstream_len) + "\n"
        return retstr

    def as_gff(self, set_name):
        """return information as a .gff line"""
        IRL = self.IRL if self.IRL != None else ""
        IRR = self.IRR if self.IRR != None else ""

        product = ""
        if self.feature and 'locus_tag' in self.feature.qualifiers:
            product = self.feature.qualifiers['locus_tag'][0]            

        annotations = ('set_id "%s"; family "%s"; group "%s"; ' + 
                       'IRL "%s"; IRR "%s"; locus_tag "%s"'
                            ) % (set_name, str(self.family),
                            str(self.group), self.IRL, self.IRR, product)
        return "\t".join([self.chromosome, "OASIS", "IS",
                str(self.start+1), str(self.end), ".",
                DIRECTION_TO_STRAND[self.direction], ".", annotations])

    def fasta(self, set_name):
        """return a fasta SeqRecord object of this IS"""
        seq_id = "|".join(["_".join((self.chromosome, str(self.start + 1),
                          str(self.end))), DIRECTION_TO_STRAND[self.direction],
                          set_name])
        if self.direction == 1:
            seq = self.seq
        else:
            seq = self.seq.reverse_complement()
        return SeqRecord.SeqRecord(id=seq_id, seq=seq, description="")

    def aa_fasta(self):
        """return a list of fasta SeqRecord objects of the ORFs of this IS"""
        seq_id = (self.chromosome + "_" +
                    str(self.start + 1) + "_" + str(self.end) + "_ORF")
        retseqs = []
        if len(self.aaseqs)==1:
            return [SeqRecord.SeqRecord(id=seq_id, seq=self.aaseqs[0],
                                description="")]
        c = 1
        for s in self.aaseqs:
            retseqs.append(SeqRecord.SeqRecord(id=seq_id+"_"+str(c), seq=s,
                                                description=""))
            c = c + 1
        #print len(retseqs)
        return retseqs
