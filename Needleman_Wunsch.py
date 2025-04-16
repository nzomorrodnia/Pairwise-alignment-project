MATCH = 1
MISMATCH = -1
GAP = -1

from score_matrix import score

def init_sides(matrix):
    for i in range(1, len(matrix)):
        matrix[i][0] = matrix[i-1][0] + GAP
    for j in range(1, len(matrix[0])):
        matrix[0][j] = matrix[0][j-1] + GAP

def rev(seq1, seq2):
    l_seq1 = list(seq1)
    l_seq2 = list(seq2)
    l_seq1.reverse()
    l_seq2.reverse()
    return ''.join(l_seq1) + '\n' + ''.join(l_seq2)

def traceback(matrix, seq1, seq2):
    aligned_seq1 = ""
    aligned_seq2 = ""

    i = len(seq1)
    j = len(seq2)
    while i > 0 or j > 0:
        if (i > 0 and j > 0) and matrix[i][j] == (matrix[i-1][j-1] + score(i,j,seq1,seq2)):
            aligned_seq1 += seq1[i-1]
            aligned_seq2 += seq2[j-1]
            i -= 1
            j -= 1
        elif i > 0 and matrix[i][j] == matrix[i-1][j] + GAP:
            aligned_seq1 += seq1[i-1]
            aligned_seq2 += "-"
            i -= 1
        else:
            aligned_seq1 += "-"
            aligned_seq2 += seq2[j-1]
            j -= 1

    return rev(aligned_seq2, aligned_seq1)

def needleman_wunsch(seq1, seq2):
    score_matrix = [[0 for _ in range(len(seq2) + 1)] for _ in range(len(seq1) + 1)]
    init_sides(score_matrix)

    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            score_matrix[i][j] = max(score_matrix[i-1][j-1] + score(i,j,seq1,seq2),score_matrix[i][j-1] + GAP, score_matrix[i-1][j] + GAP)
    for row in score_matrix:
        print(row)
    return traceback(score_matrix, seq1, seq2)


print(needleman_wunsch("AGCCA", "AGGACT"))
# AGC-CA
# ++--+- = 0
# AGGACT

# A-GCCA
# +-+-+- = 0
# AGGACT

print(needleman_wunsch("GATTACA","GCATGCG"))
# GCATG-CG
# +-++--+- = 0
# G-ATTACA

# GCA-TGCG
# +-+-+-+- = 0
# G-ATTACA

print(needleman_wunsch("ATTGC", "AGGC"))
# ATTGC
# +--++ = 1
# A-GGC

# ins1 = "ATGGCCTTCTGGCTCCAAGCTGCATCTCTGCTGGTGTTGCTGGCGCTCTCCCCCGGGGTAGATGCTGCAGC" \
# # "TGCCCAGCACCTGTGTGGCTCTCACCTGGTGGACGCCCTCTATCTGGTGTGTGGAGAGAAAGGATTCTTTTACACCCCAAAGAGAGATGT" \
# # "GGATCCCCTTATAGGGTTCCTCTCTCCAAAATCAGCAAAGGAGAACGAAGAGTACCCCTTCAAAGACCAGACGGAGATGATGGTAAAGAGAGGTATTGTAGA" \
# # "GCAGTGCTGTCACAAGCCCTGCAACATCTTCGACCTGCAAAACTACTGCAACTGA"

# ins2 = "ATGGCCCTGTGGATGCACCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCCGAGCCAGC" \
# # "CCCGGCCTTTGTGAACCAGCACCTGTGCGGCCCCCACCTGGTGGAAGCCCTCTACCTGGTGTGCGGGGAGCGAGGTT" \
# # "TCTTCTACGCACCCAAGACCCGCCGGGAGGCGGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGTGGGGGCTCTATCACGG" \
# # "GCAGCCTGCCACCCTTGGAGGGTCCCATGCAGAAGCGTGGCGTCGTGGATCAGTGCTGCACCAGCATCTGCTCCCTCTACCAGCTGCAGAACTACTGCAACTAG"

# print(needleman_wunsch(ins1,ins2))
