MATCH = 1
MISMATCH = -1
GAP = -10

from score_matrix import score

def init_sides(matrix):
    for i in range(1, len(matrix)):
        matrix[i][0] = matrix[i-1][0] + GAP
    for j in range(1, len(matrix[0])):
        matrix[j][0] = matrix[j-1][0] + GAP


def needleman_wunsch(seq1, seq2):
    score_matrix = [[0 for _ in range(len(seq2) + 1)] for _ in range(len(seq1) + 1)]
    init_sides(score_matrix)
    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):
            score_matrix[i][j] = max(score_matrix[i-1][j-1] + score(i,j,seq1,seq2), score_matrix[i-1][j] + GAP, score_matrix[i][j-1] + GAP)
    for row in score_matrix:
        print(row)

print(needleman_wunsch("GATTACA", "GCATGCG"))
