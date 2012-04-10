"""EdgeFinder.py
Created by David Robinson
6/10/09

Finds the edge of a conserved region that is repeated using the maximum
likelihood estimate, assuming a multinomial distribution
"""

### IMPORTS ###
import random
import math

from OASIS_functions import *
from Constants import *

### FIXED PARAMETERS ###
# distribution of nucleotides if there is no conservation
nucs = ["A", "G", "C", "T"]

NULL_DISTRIBUTION = dict([(n, .25) for n in nucs])

# distribution of nucleotides if there is conservation
# assume that the most common nucleotide has a large probability
# and that each other a small probability
[MOST_COMMON, OTHER] = [.991, .003]

# nucleotide distributions- map nucleotides to parameter distributions
CONSERVED_DISTRIBUTIONS = dict([(n, dict([(n, OTHER) for n in nucs])) for n in nucs])
for n in nucs:
    CONSERVED_DISTRIBUTIONS[n][n] = MOST_COMMON

REMEMBERED_FACTORIALS = {}

### FUNCTIONS ###
def mean(lst):
    return sum(lst) / float(len(lst))

def factorial(n):
    """factorial of an integer"""
    # remember factorial for future
    if n in REMEMBERED_FACTORIALS:
        return REMEMBERED_FACTORIALS[n]
    result = 1
    for i in xrange(1, abs(n)+1):
        result *= i
    if n < 0:
        result = -result
    REMEMBERED_FACTORIALS[n] = result
    return result

def product(lst):
    """return the product of all elements of a list"""
    return reduce(lambda x, y: x * y, lst, 1)

def mode_nucs(lst):
    """return the most common nucleotide of a list"""
    return argmax(nucs, lambda e: lst.count(e))

def argmin(seq, fn):
    """Return an element with lowest fn(seq[i]) score; tie goes to first one.
    >>> argmin(['one', 'to', 'three'], len)
    'to'
    """
    best = seq[0]; best_score = fn(best)
    for x in seq:
        x_score = fn(x)
        if x_score < best_score:
            best, best_score = x, x_score
    return best

def argmax(seq, fn):
    """Return an element with highest fn(seq[i]) score; tie goes to first one.
    >>> argmax(['one', 'to', 'three'], len)
    'three'
    """
    return argmin(seq, lambda x: -fn(x))

def multinomial_likelihood(data, parameters):
    """given a list of data and a dictionary mapping values to probabilities,
    return the probability of this data"""
    n = len(data)
    
    # get counts
    counts = [(v, data.count(v)) for v in parameters]
    
    term1 = factorial(n) / product([factorial(c) for v, c in counts])
    term2 = product([math.pow(parameters[v], c) for v, c in counts])
    return term1 * term2

def log_multinomial_likelihood(data, parameters):
    """given a list of data and a dictionary mapping values to probabilities,
    return the probability of this data"""
    n = len(data)
    
    # get counts
    counts = [(v, data.count(v)) for v in parameters]
    
    term1 = math.log(factorial(n) / product([factorial(c) for v, c in counts]))
    term2 = sum([math.log(math.pow(parameters[v], c)) for v, c in counts])
    return term1 + term2

def log_likelihood_conserved(data):
    """return the log likelihood that this set of nucleotides is conserved"""
    # find the mode
    most_common = mode_nucs(data)
    
    # find the likelihood given that the mode is the most common parameter
    # and all others are less common
    return log_multinomial_likelihood(data, CONSERVED_DISTRIBUTIONS[most_common])

def log_likelihood_unconserved(data):
    """return the likelihood that this set of nucleotides is from an unconserved
    region"""
    return log_multinomial_likelihood(data, NULL_DISTRIBUTION)

def find_left_edge(seqs):
    """given a list of aligned Seq objects, find the index of the edge
    of conservation among them assuming that the left side is conserved"""
    # makes sure there are multiple sequences and that sequences are the same length
    assert len(seqs) > 1, "find_left_edge must be given multiple sequences"
    
    #print "\n".join(map(str, seqs))
    
    lengths = list(set([len(s) for s in seqs]))
    
    if 0 in lengths:
        # missing region- return 0 as distance out
        return 0
    
    if len(lengths) != 1:
        print "there are", len(lengths), "possible lengths!"
        print lengths
        return 0
    #assert len(lengths) == 1, "All find_left_edge sequences must be the same length"
    
    n = lengths[0]
    
    # get distributions of nucleotides at each index
    dists = [[s[i] for s in seqs] for i in range(n)]
    
    # maximize the likelihood
    # this is the product of the likelihoods that it is conserved up to the
    # nth element and that it is not conserved from the nth element onward
    
    # compute both LOG likelihoods for each element
    conserved_l = [log_likelihood_conserved(d) for d in dists]
    unconserved_l = [log_likelihood_unconserved(d) for d in dists]
    
    # dynamic programming
    # find the likelihoods of all possible lengths of conserved sequences
    # starting from the beginning of the sequence, and all possible unconserved
    # leading to the end of the sequence
        
    conserved_region_l = [0] * (n + 1)
    log_likelihood = 0
    for i in range(n):
        conserved_region_l[i] = log_likelihood
        log_likelihood = log_likelihood + conserved_l[i]
    conserved_region_l[i] = log_likelihood
    
    # same for likelihoods of corresponding unconserved sequences (going backwards
    unconserved_region_l = [0] * (n + 1)
    log_likelihood = 0
    for i in range(n, -1, -1):
        unconserved_region_l[i] = log_likelihood
        if i < n:
            log_likelihood = log_likelihood + unconserved_l[i]
    
    #for k in range(n):
    #    print "Likelihood", k, product(conserved_l[:k]) * product(unconserved_l[k:])
    
    # find k with highest log likelihood:
    highest_lik = None
    highest_k = -1
    for k in xrange(0, n):
        lik = conserved_region_l[k] + unconserved_region_l[k]
        if lik > highest_lik or k == 0:
            highest_lik = lik
            highest_k = k
        
    return highest_k

def find_right_edge(seqs):
    """find the index edge of conservation, this time assuming the right is
    conserved and the left is unconserved"""
    backwards_seqs = [s[::-1] for s in seqs]
    return len(seqs[0]) - find_left_edge(backwards_seqs)

### CLASSES ###
class ConservedElement:
    """a suspected conserved element, with a known conserved region and
    surrounding unknown regions (each a Seq object)"""
    def __init__(self, known, leftregion, rightregion, knownstart=None, knownend=None):
        self.known = known
        self.leftregion = leftregion
        self.rightregion = rightregion
        
        # the known start and end positions of the region
        self.knownstart = knownstart
        self.knownend = knownend

class ConservedElementSet:
    """a set of elements that have some conserved region"""
    def __init__(self, lst=None):
        """can be initialized with a list of elements or with nothing"""
        self.elements = []
        
        if lst:
            self.elements = self.elements + lst[:] # copy of list
    
    def add_element(self, element):
        """add a ConservedElement"""
        self.elements.append(element)
    
    def find_edges(self):
        """find and mark the edges"""
        # use maximum likelihood edges
        self.left_side_edge = find_right_edge([e.leftregion for e in self.elements])
        self.right_side_edge = find_left_edge([e.rightregion for e in self.elements])
        
        # figure out the offset for the left and right sides
        self.left_offset = self.left_side_edge - len(self.elements[0].leftregion)
        self.right_offset = self.right_side_edge
        
        for e in self.elements:
            e.knownstart = e.knownstart + self.left_offset
            e.knownend = e.knownend + self.right_offset
    
    def __repr__(self):
        """for testing purposes"""
        ret = ""
        for e in self.elements:
            ret = ret + str(e.leftregion[:self.left_side_edge]) + "|" + \
            str(e.leftregion[self.left_side_edge:]) + " " + str(e.known) + " " + \
            str(e.rightregion[:self.right_side_edge]) + "|" + \
            str(e.rightregion[self.right_side_edge:]) + " " + str(e.knownstart) + ", " + str(e.knownend) + "\n"
        return ret

### MAIN ###

# experiment functions
def random_seq(length):
    return "".join([random.choice(nucs) for i in range(length)])

# experiment
if __name__ == "__main__" and False:
    real_element = "AAAAAGAATACGGGGAAAAA"
    element_regions = [(random_seq(10) + real_element[:10], real_element[4:-6], real_element[-15:] + random_seq(10)) for i in range(4)]
    elements = [ConservedElement(k, l, r, 223, 233) for l, k, r in element_regions]
    
    s = ConservedElementSet(elements)
    s.find_edges()
    print s