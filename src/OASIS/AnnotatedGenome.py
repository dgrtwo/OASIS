"""AnnotatedGenome.py
Created by David Robinson
8/11/08

OASIS is a module designed for IS element annotation for prokaryote
genomes.

This module represents the data on a genome that allows it to be annotated.
This data can be imported from a genbank file or given to it as a list of
SeqRecords."""

#imports
import os
import glob

from Bio import Seq
from Bio import SeqIO
from Bio.Blast import NCBIStandalone
from Bio.Blast import NCBIXML

import ISSet
import IS
import Profile
from OASIS_functions import *
from Constants import *

#classes
class AnnotatedGenome:
    """The genome class represents a genome and data on its gene that is needed
    to annotate the IS elements of a genome"""

    def __init__(self, genome_file=None, lst=None, annotated=False,
                 data_folder=None):
        """initializes with an genbank seq_record object and whether
        the genome's transposases have been annotated"""
        #transform the genome into a dictionary
        if genome_file:
            inf = open(genome_file)
            seqlist = list(SeqIO.parse(inf, "genbank"))
            inf.close()

            self.genome_file = genome_file
            #seqlist = list(genbankparse(genome_file))

            # REMEMBER SEQUENCES
            self.seqlist = seqlist
        elif lst:
            seqlist = lst
        else:
            return False #no data given

        self.namedict = {}
        self.seqdict = SeqIO.to_dict(seqlist)
        self.candidates = []
        self.annotations = None

        #set up an ISFinder profile
        #self.profile = Profile.Profile("corrected_linker_out.fasta", "corrected_tpase_out.fasta")
        if data_folder == None:
            data_folder = os.path.join(os.path.split(__file__)[0], "data")
        aafile = os.path.join(data_folder, ISFINDER_AA_FILE)
        self.profile = Profile.Profile(aafile)

        for seq_record in seqlist:
            if annotated:
                #find all transposase genes
                thislist = [f for f in seq_record.features if related(f) and 'locus_tag' in f.qualifiers]
                for f in thislist:
                        self.namedict[f.qualifiers['locus_tag'][0]] = seq_record.id
                self.candidates = self.candidates + thislist

                #if any([related(f) for f in seq_record.features]):
                #    print seq_record
            else:
                #include all genes
                thislist = [f for f in seq_record.features if f.type == "CDS"]
                for f in thislist:
                    self.namedict[f.qualifiers['locus_tag'][0]] = seq_record.id
                self.candidates = self.candidates + thislist

    def annotate(self):
        """fill self.annotations with IS elements in this genome"""
        #if it's already been annotated, stop
        if self.annotations: return

        self.annotations = []
        already_indices = []
        already = []

        list_length = len(self.candidates)

        #iterate through candidates, compare each to the rest in the list
        for i in range(list_length-1):
            amatch = None
            lst = []
            for j in range(i+1, list_length):
                if j not in already_indices: #no need to compare more than once
                    if self.__match(i, j):
                        lst.append(self.__find_full(i, j))

                        # check if there are multiple different lengths
#                         length_i = self.candidates[i].location._end.position - self.candidates[i].location._start.position
#                         length_j = self.candidates[j].location._end.position - self.candidates[j].location._start.position
#                         if length_i != length_j:
#                             print self.candidates[i]
#                             print self.candidates[j]

                        amatch = self.__find_full(j, i)
                        already_indices.append(j)
                        if i not in already_indices:
                            already_indices.append(i)

            #remove multiple ORFs
            if amatch:
                lst = [amatch] + lst
                newlst = [m for m in lst if not IS_in_list(m, already)]
                already = already + newlst
                if len(newlst) > 0:
                    self.annotations.append(ISSet.ISSet(newlst, self.profile))

        #perform operations to improve IS elements
        self.clean_up()

    def clean_up(self):
        """perform the necessary changes- clean edges, add IRs, filter,
        BLAST, find partials and singles, etc."""
        self.annotations = [is_set for is_set in self.annotations if is_set.filter()]
        #[is_set.clean_edges() for is_set in self.annotations]
        [is_set.re_annotate() for is_set in self.annotations]
        #print "a"
        self.__find_singles()
        self.__find_partials(minimum_blast_length=500)
        #print "b"
        [is_set.add_IRs() for is_set in self.annotations]
        [is_set.identify_transposase() for is_set in self.annotations]
        #print "c"
        self.__group_annotations()
        [is_set.re_annotate(from_edges=True) for is_set in self.annotations]
        return

        #below here doesn't happen
        self.__find_partials()
        self.__find_singles()
        [is_set.clean_edges() for is_set in self.annotations]
        [is_set.add_IRs() for is_set in self.annotations]
        self.__find_partials()
        [is_set.clean_edges() for is_set in self.annotations]
        [is_set.add_IRs() for is_set in self.annotations]
        self.remove_redundant()

    def __group_annotations(self):
        """group self.annotations by combining sets of similar IS elements"""
        i = 0
        n = len(self.annotations)
        marked = []
        while i < n-1:
            j = i + 1
            while j < n:
                if similar(self.annotations[i], self.annotations[j]):
                    #print "i =", i, "j =", j,
                    if self.annotations[i].score() > self.annotations[j].score():
                    #    print "adding j to marked"
                        marked.append(j)
                    else:
                    #    print "adding i to marked"
                        marked.append(i)
                j = j + 1
            i = i + 1

        #print "MARKED", marked

        #for i, a in enumerate(self.annotations):
        #    print "Annotation", i
        #    print a

        self.annotations = [a for a, i in zip(self.annotations, range(n)) if i not in marked]

    def write_annotations(self, outfilename, folder=None):
        """given an output file, write the annotations out to a text file and
        a fasta file"""
        #get rid of the outfile's extension and any slashes
        outfile = os.path.splitext(outfilename)[0]


        #try:
        #if there is a folder, add it
        if folder: outfile = os.path.join(folder, outfile)

        #annotate the genome
        self.annotate()

        #write files
        self.write_txt(outfile)
        self.write_fasta(outfile)

        # print confirmation
        print "Completed", outfile

    def write_txt(self, outfile):
        """write annotations as a gff file"""
        outf = open(outfile + ".gff", "w")
        for i, is_set in enumerate(self.annotations):
            outf.write("\n".join(is_set.as_gff(i + 1)) + "\n")
        outf.close()

    def write_fasta(self, outfile):
        """write fasta file"""
        fasta_out = open(outfile + ".fasta", "w")
        for is_set in self.annotations:
            SeqIO.write(is_set.fasta(), fasta_out, "fasta")
            fasta_out.write("\n")
        fasta_out.close()

    def __write_singles(self, outfile):
        """write a single IS from each group to the output file"""
        fasta_out = open(outfile, "w")
        for is_set in self.annotations:
            SeqIO.write([is_set.representative().fasta()], fasta_out, "fasta")
        fasta_out.close()

    def __get_extensions(self, index, window):
        """given the index of a candidate, get the windows surrounding it"""
        f = self.candidates[index]
        start = f.location._start.position
        end = f.location._end.position
        chromosome = self.namedict[f.qualifiers['locus_tag'][0]]
        seq = self.seqdict[chromosome].seq
        if self.candidates[index].strand == None: #if it's forward
            return seq[start-window:start].reverse_complement(), seq[end:end+window]
        return seq[end:end+window], seq[start-window:start].reverse_complement()

    def __match(self, i, j):
        """test whether two CDSs match as IS gene candidates"""
        #first, check whether their CDSs are the same length and have low
        #divergence

        aaseq_1 = self.get_aaseq(self.candidates[i])
        aaseq_2 = self.get_aaseq(self.candidates[j])

        if abs(len(aaseq_1) - len(aaseq_2)) > LENGTH_WIGGLE:
            return False #they are not the same length

        if divergence(str(aaseq_1), str(aaseq_2)) < MIN_DIVERGENCE:
            return False #the amino acid sequences are too divergent

        #now, check their surrounding regions
        #get windows
        before_window_1, after_window_1 = self.__get_extensions(i, EXTEND_LENGTH)
        before_window_2, after_window_2 = self.__get_extensions(j, EXTEND_LENGTH)

        #check whether they match the similarity threshold
        if not seq_match(before_window_1, before_window_2, EXTEND_LENGTH, MAX_MISS) \
            and not seq_match(after_window_1, after_window_2, EXTEND_LENGTH, MAX_MISS):
            return False
        #if they have gotten this far, they are a good match
        return True

    def __find_full(self, i, j):
        """given the indices of two genes, return a match object for the
        IS containing the second gene that matches the first"""
        feature = self.candidates[j]
        before_window_1, after_window_1 = self.__get_extensions(i, FULL_EXTEND)
        before_window_2, after_window_2 = self.__get_extensions(j, FULL_EXTEND)

        start = feature.location._start.position
        end = feature.location._end.position

        #correct for extension length
        if feature.strand == None: #if it is direct
            start = start - extend(before_window_1, before_window_2, FULL_MAX_MISS)
            end = end + extend(after_window_1, after_window_2, FULL_MAX_MISS)
        else: #if it is reversed
            start = start - extend(after_window_1, after_window_2, FULL_MAX_MISS)
            end = end + extend(before_window_1, before_window_2, FULL_MAX_MISS)
        chromosome = self.namedict[feature.qualifiers['locus_tag'][0]]
        return IS.IS(feature, chromosome, start, end, self)

    def __find_partials(self, minimum_blast_length=0):
        """find partial IS elements by blasting the sequences against the
        genome"""
        #if there are no IS elements, skip this step
        if len(self.annotations) == 0: return

        #write a temporary genome fasta file
        blast_db = os.path.join(TEMPORARY_DIRECTORY, "OASIS_temp_genome.fasta")
        outf = open(blast_db, "w")
        SeqIO.write(self.as_records(), outf, "fasta")
        outf.close()
        #turn it into a database
        os.system(FORMAT_EXE + " -p F -i " + blast_db)

        #write a temporary IS fasta file
        blast_file = os.path.join(TEMPORARY_DIRECTORY, "OASIS_temp_IS.fasta")
        self.__write_singles(blast_file)

        #get the directions of these sample IS's
        directions = [is_set.lst[0].direction for is_set in self.annotations]

        #clear annotations
        self.annotations = []

        #perform a blast
        result_handle, error_handle = NCBIStandalone.blastall(BLAST_EXE,
                                        "blastn", blast_db, blast_file)
        blast_records = NCBIXML.parse(result_handle)

        #iterate over the results and the directions of the queries
        for record, sample_direction in zip(blast_records, directions):
            ISlist = []
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    if hsp.expect < E_VALUE_CUTOFF and len(hsp.sbjct) >= MIN_PARTIAL_LEN and \
                        len(hsp.sbjct) > minimum_blast_length:
                        chromosome = alignment.title.split(" ")[1]
                        start = hsp.sbjct_start-1
                        end = start + len(hsp.sbjct)
                        #find out what the gene is
                        f = self.get_feature(chromosome, start, end)
                        thisdir = hsp.frame[1] * sample_direction
                        ISlist.append(IS.IS(f, chromosome, start, end, self, dir=thisdir))
            if len(ISlist) > 0:
                self.annotations.append(ISSet.ISSet(ISlist, self.profile))

        #clean up- remove the temporary files
        os.remove(blast_db)
        os.remove(blast_file)
        for f in glob.glob(blast_db + ".n*"):
            os.remove(f)
        os.remove("formatdb.log")

    def get_aaseq(self, f):
        """given a feature, get the amino acid sequence"""
        # if possible, use given translation
        if "translation" in f.qualifiers:
            return Seq.Seq(f.qualifiers['translation'][0])

        # otherwise, translate from genome
        id = f.qualifiers['locus_tag'][0]
        start = f.location._start.position
        end = f.location._end.position
        thisseq = self.seqdict[self.namedict[id]].seq[start:end]
        if f.strand == -1:
            thisseq = thisseq.reverse_complement()
        return translate(thisseq)

    def get_chrom(self, chromosome):
        """given the name of a chromosome, return its sequence"""
        return self.seqdict[chromosome].seq

    def as_records(self):
        """return this genome as a list of SeqRecord objects"""
        return self.seqdict.values()

    def get_feature(self, chromosome, start, end):
        """given a chromosome, a start, and an end, return a feature in that
        segment, or None if not"""
        for f in self.candidates:
            f_chrom = self.namedict[f.qualifiers['locus_tag'][0]]
            if chromosome != f_chrom:
                continue
            f_start = f.location._start.position
            f_end = f.location._end.position
            if f_start > start-FEATURE_WIGGLE and f_end < end+FEATURE_WIGGLE:
                return f

    def get_multiple_features(self, chromosome, start, end):
        f_lst = []
        for f in self.candidates:
            f_chrom = self.namedict[f.qualifiers['locus_tag'][0]]
            if chromosome != f_chrom:
                continue
            f_start = f.location._start.position
            f_end = f.location._end.position
            if f_start > start-FEATURE_WIGGLE and f_end < end+FEATURE_WIGGLE:
                f_lst.append(f)

        return f_lst

    def __find_singles(self):
        """add single copy IS elements to the annotations"""
        #iterate over IS elements that haven't already been used (no multi ISs)
        unused = self.__unused()

        i = 0
        #iterate over unused candidates
        while i < len(unused):
            ORF_list = None
            c = unused[i]
            #check if the next is the same gene
            if i < len(unused) - 1 and same_gene(c, unused[i+1]):
                ORF_list = [c, unused[i+1]]
                i = i + 1
            if not self.redundant_feature(c):
                result_IS = None
                if ORF_list:
                    result_IS = self.find_IRs_around(lst=ORF_list)
                else:
                    result_IS = self.find_IRs_around(c=c)
                if result_IS:
                    self.annotations.append(ISSet.ISSet([result_IS], self.profile))
            i = i + 1

    def all_elements(self):
        """return a list of all IS elements"""
        ret = []
        for set in self.annotations:
            for e in set.lst:
                ret.append(e)
        return ret

    def __unused(self):
        """return a list of candidates that have not yet been annotated"""
        return [c for c in self.candidates if not self.redundant_feature(c)]

    def redundant_feature(self, f):
        """check whether a feature is already contained within an annotated IS
        element"""
        chromname = self.namedict[f.qualifiers['locus_tag'][0]]
        start = f.location._start.position
        end = f.location._end.position
        for e in self.all_elements():
            if e.chromosome == chromname and \
                e.start - REDUNDANT_WIGGLE < start and \
                e.end + REDUNDANT_WIGGLE > end:
                    return True
        return False

    def redundant_IS(self, e):
        """check whether an IS element is redundant"""
        return IS_in_list(e, self.all_elements())

    def find_IRs_around(self, c=None, lst=None):
        """given a candidate gene, or list of ORFs, return an IS around it,
        or return False if no inverted repeats are found"""
        if c:
            seq = self.get_aaseq(c)
            start = c.location._start.position
            end = c.location._end.position
            dir = {None: 1, 1: 1, -1:-1}[c.strand]
            chromname = self.namedict[c.qualifiers['locus_tag'][0]]
        elif lst:
            seq = self.get_aaseq(lst[0])
            start = min([c.location._start.position for c in lst])
            end = max([c.location._end.position for c in lst])
            dir = {None: 1, 1: 1, -1:-1}[lst[0].strand]
            chromname = self.namedict[lst[0].qualifiers['locus_tag'][0]]
            c = lst[0]
        else:
            return False

        chromosome = self.get_chrom(chromname)

        if dir == 1:
            before = chromosome[start-SINGLE_IR_OUT_WINDOW:start+SINGLE_IR_IN_WINDOW]
            after = chromosome[end-SINGLE_IR_IN_WINDOW:end+SINGLE_IR_OUT_WINDOW]
            after = after.reverse_complement()
        else:
            before = chromosome[end-SINGLE_IR_IN_WINDOW:end+SINGLE_IR_OUT_WINDOW]
            after = chromosome[start-SINGLE_IR_OUT_WINDOW:start+SINGLE_IR_IN_WINDOW]
            before = before.reverse_complement()
        #IR_result = find_IR_large_window(before, after.reverse_complement())

        #family, group = self.profile.identify_family(seq)
        IR_result = self.profile.find_IRs(None, before, after, SINGLE_IR_IN_WINDOW)

        if IR_result:
            IR1, IR2 = IR_result
            full_window = chromosome[start-SINGLE_IR_OUT_WINDOW:end+SINGLE_IR_OUT_WINDOW]

            try:
                if dir == 1:
                    pos1 = str(full_window).index(str(IR1))
                    pos2 = str(full_window).index(str(IR2))
                else:
                    pos1 = str(full_window).index(str(IR1.reverse_complement()))
                    pos2 = str(full_window).index(str(IR2.reverse_complement()))
            except ValueError:
                return False

            newstart = start - SINGLE_IR_OUT_WINDOW + pos1
            newend = start - SINGLE_IR_OUT_WINDOW + pos2 + len(IR2)
            newIS = IS.IS(c, chromname, newstart, newend, self)
            if newIS.length < MINLENGTH: return False
            return newIS

            #old stuff:
            start_index, end_index = IR_result
            print start_index, end_index
            #print start_index, end_index
            if dir == 1:
                newstart = start - start_index
                newend = end + end_index
            else:
                print after[end_index:end_index+15]
                newstart = start - end_index
                newend = end + start_index
            newIS = IS.IS(c, chromname, newstart, newend, self)
            if newIS.length < MINLENGTH: return False
            print newIS.seq
            raise Exception("Bad Sequence")
            return newIS

        return False

    def remove_redundant(self):
        """remove redundant IS annotations"""
        old_annotations = self.annotations
        self.annotations = []
        for set in old_annotations:
            if set.remove_redundant(self):
                self.annotations.append(set)
