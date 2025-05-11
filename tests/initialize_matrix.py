import numpy as np

def initialize_matrix(seq1, seq2, gap_penalty=-1):

    m = len(seq1) # the number of rows + 1
    n = len(seq2) # the number of columns + 1

    # Creating an m+1 * n+1 matrix with zeros
    matrix = np.zeros((m+1, n+1), dtype=int)

    # Filling the first row (gap penalties)
    for i in range(1, n+1): # starting range from 1 to number of columns+1
        matrix[0][i] = matrix[0][i-1] + gap_penalty

    
    # Filling the first column (gap penalties)
    for j in range(1, m+1):
        matrix[j][0] = matrix[j-1][0] + gap_penalty
    
    return matrix

# Example usage
seq1 = "GATTACA"
seq2 = "GCATGCU"
matrix = initialize_matrix(seq1, seq2)
print(matrix)