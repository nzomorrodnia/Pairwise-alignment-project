import numpy as np

MATCH = 2
MISMATCH = -1
GAP_OPEN_PENALTY = -3
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
    Q = np.zeros((n + 1, m + 1))
    P = np.zeros((n + 1, m + 1))

    for i in range(1, n + 1):
        Q[i][0] = float('-INF')
        score_matrix[i][0] = GAP_OPEN_PENALTY + i*GAP_EXTEND_PENALTY

    for j in range(1, m + 1):
        P[0][j] = float('-INF')
        score_matrix[0][j] = GAP_OPEN_PENALTY + j*GAP_EXTEND_PENALTY

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score_match = score_matrix[i - 1][j - 1] + (MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH)
            P[i][j] = max(P[i - 1][j] + GAP_EXTEND_PENALTY, score_matrix[i - 1][j] + GAP_OPEN_PENALTY + GAP_EXTEND_PENALTY)
            Q[i][j] = max(Q[i][j - 1] + GAP_EXTEND_PENALTY, score_matrix[i][j - 1] + GAP_OPEN_PENALTY + GAP_EXTEND_PENALTY)
            score_matrix[i][j] = max(score_match, P[i][j], Q[i][j])
    for row in score_matrix:
        print(row)
    for row in P:
        print(row)
    for row in Q:
        print(row)
    return traceback(seq1, seq2, score_matrix, Q, P)


def alignment_with_prints(seq1, seq2):
    alignments = needleman_wunsch(seq1, seq2)
    symbols = ""
    for n in range(len(alignments[0])):
        if alignments[0][n] == alignments[1][n]:
            symbols += "|"
        elif alignments[0][n] == '-' or alignments[1][n] == '-':
            symbols += ' '
        else:
            symbols += '*' 
    print(f"for: '{seq1} + ' and: ' + {seq2}'")
    print(f"got:\n{alignments[0]}\n{symbols}\n{alignments[1]}")
    return alignments


alignment_with_prints("CCGA","CG")
