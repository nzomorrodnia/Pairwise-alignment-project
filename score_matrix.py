import blosum as bl
import pandas as pd


def score(i,j,seq1,seq2):
    '''Calculates the BLOSUM62 score for substitution in i,j'''    
    blosum = bl.BLOSUM(62)
    key = seq1[i-1]
    for inner_key in blosum[key]:
        if inner_key  == seq2[j-1]:
            return blosum[key][inner_key]

print(bl.BLOSUM(90)['C']['D'])