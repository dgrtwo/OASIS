"""EdgeFinder.py
Created by David Robinson
6/10/09

Finds the edge of a conserved region that is repeated using the maximum
likelihood estimate, assuming a multinomial distribution
"""

### IMPORTS ###
import random
import math

from Bio import Seq

from OASIS_functions import *
from Constants import *

### FIXED PARAMETERS ###
# distribution of nucleotides if there is no conservation
nucs = ["A", "G", "C", "T"]

NULL_DISTRIBUTION = dict([(n, .25) for n in nucs])

# distribution of nucleotides if there is conservation
# assume that the most common nucleotide has a large probability
# and that each other hqs a small probability
[MOST_COMMON, OTHER] = [.995, .005]

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

def cumulative_sum(lst):
    """
    return the cumulative sum of a list (will end up being n+1 long and start
    with 0)
    """
    ret = [0]
    total = 0
    for e in lst:
        total += e
        ret.append(total)
    return ret

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

def whichmax(seq):
    """return the index of the highest element"""
    return max([(b, a) for a, b in enumerate(seq)])[1]

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

def combine_likelihoods(conserved_likelihoods, unconserved_likelihoods):
    conserved_region_l = cumulative_sum(conserved_likelihoods)
    unconserved_region_l = cumulative_sum(unconserved_likelihoods[::-1])[::-1]
    return map(sum, zip(conserved_region_l, unconserved_region_l))

def find_left_edges(seqs, individual=True, verbose=False):
    """given a list of aligned Seq objects, find the index of the edge
    of conservation among them assuming that the left side is conserved"""
    # makes sure there are multiple sequences and that sequences are the same length
    assert len(seqs) > 1, "find_left_edge must be given multiple sequences"

    lengths = list(set([len(s) for s in seqs]))

    if 0 in lengths:
        # missing region- return 0 as distance out
        return [0] * len(seqs)

    n = lengths[0]

    # get distributions of nucleotides at each index
    dists = [[s[i] for s in seqs] for i in range(n)]

    # maximize the likelihood
    # this is the product of the likelihoods that it is conserved up to the
    # nth element and that it is not conserved from the nth element onward

    conserved_l = map(log_likelihood_conserved, dists)
    unconserved_l = map(log_likelihood_unconserved, dists)

    likelihoods = combine_likelihoods(conserved_l, unconserved_l)
    best_k = whichmax(likelihoods)

    revcmp = lambda s: str(Seq.Seq(s).reverse_complement())

    if verbose:
        print "\n".join([revcmp(s[:best_k+100]) for s in seqs])
        print "\n".join([" " * 100 + revcmp(s[:best_k]) for s in seqs])
        print best_k

    if len(seqs) < 3 or not individual:
        return [best_k] * len(seqs)

    modes = map(mode_nucs, dists[:best_k])

    individual_limits = []
    for s in seqs:
        conserved_l = [math.log(MOST_COMMON) if n == mode else math.log(OTHER)
                            for n, mode in zip(s, modes)]
        unconserved_l = [math.log(NULL_DISTRIBUTION[n]) for n in s[:best_k]]
        likelihoods = combine_likelihoods(conserved_l, unconserved_l)

        individual_limits.append(whichmax(likelihoods))

    # if there is a maximum, cut it to the next highest (it's an artifact)
    individual_limits[whichmax(individual_limits)] = (
                                        sorted(individual_limits)[-2])

    return individual_limits

def find_right_edges(seqs, individual=True, verbose=False):
    """find the index edge of conservation, this time assuming the right is
    conserved and the left is unconserved"""
    backwards_seqs = [s[::-1] for s in seqs]
    return [len(seqs[0]) - e
                for e in find_left_edges(backwards_seqs, individual, verbose)]

