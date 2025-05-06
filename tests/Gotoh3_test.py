import numpy as np


MATCH = 2
MISMATCH = -1
GAP_OPEN = -3
GAP_EXTEND = -1


def get_matrices(n, m):
    S = np.zeros((n + 1, m + 1))
    P = np.zeros((n + 1, m + 1))
    Q = np.zeros((n + 1, m + 1))

    S[0][0] = 0
    P[0][0] = float('NaN')
    Q[0][0] = float('NaN')

    for i in range(1, n+1):
        S[i][0] = GAP_OPEN + (i) * GAP_EXTEND if i == 1 else S[i - 1][0] + GAP_EXTEND
        Q[i][0] = float('-inf')
    for j in range(1, m+1):
        S[0][j] = GAP_OPEN + (j) * GAP_EXTEND if j == 1 else S[0][j - 1] + GAP_EXTEND
        P[0][j] = float('-inf')

    return S, P, Q


def traceback(seq1, seq2, S, P, Q, i, j, current_alignment, alignments, matrix_to_use):
    if i == 0 and j == 0:
        alignments.append(current_alignment)
        return

    '''Issue when only one of either i or j is 0 - try using != instead'''
    if matrix_to_use == 'S':
        score = S[i][j]
        if i > 0 and j > 0 and S[i - 1][j - 1] + (MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH) == score:
            traceback(seq1, seq2, S, P, Q, i - 1, j - 1, (current_alignment[0] + seq1[i - 1], current_alignment[1] + seq2[j - 1]), alignments, 'S')
        if i > 0 and j > 0 and P[i-1][j-1] + (MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH) == score:
            traceback(seq1, seq2, S, P, Q, i - 1, j - 1, (current_alignment[0] + seq1[i - 1], current_alignment[1] + seq2[j - 1]), alignments, 'P')
        if i > 0 and j > 0 and Q[i-1][j-1] + (MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH) == score:
            traceback(seq1, seq2, S, P, Q, i - 1, j - 1, (current_alignment[0] + seq1[i - 1], current_alignment[1] + seq2[j - 1]), alignments, 'Q')
        if i > 0 and P[i][j] == score:
            traceback(seq1, seq2, S, P, Q, i - 1, j, (current_alignment[0] + seq1[i - 1], current_alignment[1] + '-'), alignments, 'P')
        if j > 0 and Q[i][j] == score:
            traceback(seq1, seq2, S, P, Q, i, j - 1, (current_alignment[0] + '-', current_alignment[1] + seq2[j - 1]), alignments, 'Q')
        if i == 0 and j > 0:
            traceback(seq1, seq2, S, P, Q, i, j - 1, (current_alignment[0] + '-', current_alignment[1] + seq2[j - 1]), alignments, 'S')
        if j == 0 and i > 0:
            traceback(seq1, seq2, S, P, Q, i - 1, j, (current_alignment[0] + seq1[i - 1], current_alignment[1] + '-'), alignments, 'S')

    elif matrix_to_use == 'P':
        score = P[i][j]
        if i > 0 and S[i - 1][j] + GAP_OPEN + GAP_EXTEND == score:
            traceback(seq1, seq2, S, P, Q, i - 1, j, (current_alignment[0] + seq1[i - 1], current_alignment[1] + '-'), alignments, 'S')
        if i > 0 and P[i-1][j] + GAP_EXTEND == score:
            traceback(seq1, seq2, S, P, Q, i - 1, j, (current_alignment[0] + seq1[i - 1], current_alignment[1] + '-'), alignments, 'P')

    elif matrix_to_use == 'Q':
        score = Q[i][j]
        if j > 0 and S[i][j - 1] + GAP_OPEN + GAP_EXTEND == score:
            traceback(seq1, seq2, S, P, Q, i, j - 1, (current_alignment[0] + '-', current_alignment[1] + seq2[j - 1]), alignments, 'S')
        if j > 0 and Q[i][j-1] + GAP_EXTEND == score:
            traceback(seq1, seq2, S, P, Q, i, j - 1, (current_alignment[0] + '-', current_alignment[1] + seq2[j - 1]), alignments, 'Q')


def all_reversed_and_no_dupes(alignments):
    alignments_rev = []
    for alignment in alignments:
        l_seq1 = list(alignment[0])
        l_seq1.reverse()
        l_seq2 = list(alignment[1])
        l_seq2.reverse()
        alignments_rev.append((''.join(l_seq1), ''.join(l_seq2)))

    return list(dict.fromkeys(alignments_rev))


def needleman_wunsch(seq1, seq2):
    n, m = len(seq1), len(seq2)
    S, P, Q = get_matrices(n, m)

    for i in range(1, n+1):
        for j in range(1, m+1):
            match_score = MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH
            P[i][j] = max(S[i-1][j] + GAP_OPEN + GAP_EXTEND, P[i - 1][j] + GAP_EXTEND)
            Q[i][j] = max(S[i][j-1] + GAP_OPEN + GAP_EXTEND, Q[i][j - 1] + GAP_EXTEND)
            S[i][j] = max(S[i-1][j-1] + match_score, P[i][j], Q[i][j])


    max_score = max(S[n][m], P[n][m], Q[n][m])
    alignments = []
    if S[n][m] == max_score:
        traceback(seq1, seq2, S, P, Q, n, m, ('', ''), alignments, 'S')
    if P[n][m] == max_score:
        traceback(seq1, seq2, S, P, Q, n, m, ('', ''), alignments, 'P')
    if Q[n][m] == max_score:
        traceback(seq1, seq2, S, P, Q, n, m, ('', ''), alignments, 'Q')
    
    for row in S:
        print(row)
    for row in P:
        print(row)
    for row in Q:
        print(row)

    return all_reversed_and_no_dupes(alignments)


def alignment_with_prints(seq1, seq2):
    '''Prints alignments from previous outputs'''
    alignments = needleman_wunsch(seq1, seq2)
    print("for: \n    '" + seq1 + "' and: '" + seq2 + "'")
    alignments_as_str = ""
    for alignment in alignments:
        alignments_as_str = alignments_as_str + "    as first sequence: '" + str(alignment[0]) + "', as second sequence: '" + str(alignment[1]) + "'\n"
    print("got: \n" + alignments_as_str)
    return alignments


alignment_with_prints("CG", "CCGA")
# alignment_with_prints("GATTACA", "CATGCATCGAC")
