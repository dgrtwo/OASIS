"""my_SW.py
Created by David Robinson
1/29/08"""

import sys
import random

def matchScore(l1, l2):
    if l1 == l2:
        return 1
    return -2


gapscore = -5

bases = ["A", "G", "C", "T"]

def __random_sequence(seqlen):
    return "".join([random.choice(bases) for i in range(seqlen)])


def __matrix_max(m):
    """given a matrix, return a tuple with the indices of the maximum element and that value"""
    max_e = -1
    max_i = 0
    max_j = 0
    for i in range(len(m)):
        for j in range(len(m[0])):
            if m[i][j] > max_e:
                max_e = m[i][j]
                max_i = i
                max_j = j
    return max_i, max_j, max_e


def __matrix_walkback(D, i, j):
    """walk back to find the beginning, return its indices"""
    max_space = max([D[i-1][j], D[i][j-1], D[i-1][j-1]])
    if D[i-1][j-1] == 0: return i, j
    direction = [D[i-1][j], D[i][j-1], D[i-1][j-1]].index(max_space)
    if direction == 0: i = i - 1
    if direction == 1: j = j - 1
    if direction == 2: i, j = i - 1, j - 1
    return __matrix_walkback(D, i, j)


def align(s1, s2):
    m = len(s1)
    n = len(s2)

    if m == 0 or n == 0:
        return 0, 0, 0, 0, 0

    D = []
    for i in range(m):
        D.append([0] * n)

    #construct matrix D
    for i in range(1, m):
        for j in range(1, n):
            D[i][j] = max([0, D[i-1][j] + gapscore, D[i][j-1] + gapscore, D[i-1][j-1] + matchScore (s1[i-1], s2[j-1])])

    max_i, max_j, score = __matrix_max(D)
    start_i, start_j = __matrix_walkback(D, max_i, max_j)
    return start_i, max_i, start_j, max_j, score
