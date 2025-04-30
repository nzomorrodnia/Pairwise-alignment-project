import numpy as np

MATCH = 2
MISMATCH = -1
GAP_OPEN_PENALTY = -2
GAP_EXTEND_PENALTY = -1


def reverse_sequences(seq1, seq2):
    l_seq1 = list(seq1)
    l_seq2 = list(seq2)
    l_seq1.reverse()
    l_seq2.reverse()
    return ''.join(l_seq1), ''.join(l_seq2)


def traceback(seq1, seq2, score_matrix, seq1_gaps_matrix, seq2_gaps_matrix):
    aligned_seq1, aligned_seq2 = "", ""
    i, j = len(seq1), len(seq2)

    while i > 0 or j > 0:
        if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i - 1][j - 1] + (MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH):
            aligned_seq1 += seq1[i - 1]
            aligned_seq2 += seq2[j - 1]
            i -= 1
            j -= 1

        elif i > 0 and score_matrix[i][j] == seq1_gaps_matrix[i][j]:
            aligned_seq1 += seq1[i - 1]
            aligned_seq2 += "-"
            i -= 1

        else:
            aligned_seq1 += "-"
            aligned_seq2 += seq2[j - 1]
            j -= 1

    return reverse_sequences(aligned_seq1, aligned_seq2)


def needleman_wunsch(seq1, seq2):
    n = len(seq1)
    m = len(seq2)

    score_matrix = np.zeros((n + 1, m + 1))
    seq1_gaps_matrix = np.zeros((n + 1, m + 1))
    seq2_gaps_matrix = np.zeros((n + 1, m + 1))

    for i in range(1, n + 1):
        seq1_gaps_matrix[i][0] = float('-INF')
        score_matrix[i][0] = seq1_gaps_matrix[i][0]

    for j in range(1, m + 1):
        seq2_gaps_matrix[0][j] = float('-INF')
        score_matrix[0][j] = seq2_gaps_matrix[0][j]

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score_match = score_matrix[i - 1][j - 1] + (MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH)
            seq1_gaps_matrix[i][j] = max(seq1_gaps_matrix[i - 1][j] + GAP_EXTEND_PENALTY, score_matrix[i - 1][j] + GAP_OPEN_PENALTY)
            seq2_gaps_matrix[i][j] = max(seq2_gaps_matrix[i][j - 1] + GAP_EXTEND_PENALTY, score_matrix[i][j - 1] + GAP_OPEN_PENALTY)
            score_matrix[i][j] = max(score_match, seq1_gaps_matrix[i][j], seq2_gaps_matrix[i][j])

    return traceback(seq1, seq2, score_matrix, seq1_gaps_matrix, seq2_gaps_matrix)


def alignment_with_prints(seq1, seq2):
    alignments = needleman_wunsch("GCATGCCAT", "CATGCATCGAC")
    print("for: '" + seq1 + "' and: '" + seq2 + "'")
    print("got: '" + alignments[0] + "' and: '" + alignments[1] + "'")
    return alignments


alignment_with_prints("GCATGCCAT", "CATGCATCGAC")
