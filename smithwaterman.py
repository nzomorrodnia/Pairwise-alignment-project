
def read_fasta(file):
    '''Reads a FASTA file and returns the sequences as a strings.'''
    sequences = [] # empty list to store sequences
    with open(file, 'r') as file:
        seq = "" # empty string
        for line in file: # read line by line
            if line.startswith('>'): # if it's a header line
                if seq: # if there's a sequence already
                    sequences.append(seq) # append to list
            else: # if not a header line
                seq += line.strip().upper() # remove whitespace and convert to uppercase
        if seq: # if there's a sequence left at the end
            sequences.append(seq) # append the last sequence
    if len(sequences) != 2:
        raise ValueError("FASTA file must contain exactly two sequences.")
    return sequences[0], sequences[1] # return the two sequences





MATCH = 1
MISMATCH = -1
GAP = -1

from score_matrix import score

def init_zero_matrix(rows, cols):
    """Initializes a score matrix filled with zeros for local alignment."""
    return [[0 for _ in range(cols)] for _ in range(rows)]

def rev(seq1, seq2):
    l_seq1 = list(seq1)
    l_seq2 = list(seq2)
    l_seq1.reverse()
    l_seq2.reverse()
    return ''.join(l_seq1) + '\n' + ''.join(l_seq2)

def traceback_sw(matrix, seq1, seq2):
    """Traceback for Smith-Waterman: starts at max value, ends at 0."""
    # Find max value in matrix
    max_score = 0
    max_pos = (0, 0)
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_pos = (i, j)

    aligned_seq1 = ""
    aligned_seq2 = ""
    i, j = max_pos

    # Traceback until score is 0
    while matrix[i][j] != 0:
        if (i > 0 and j > 0 and 
            matrix[i][j] == matrix[i-1][j-1] + score(i, j, seq1, seq2)):
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

    return rev(aligned_seq2, aligned_seq1)  # reverse & return nicely

def smith_waterman(seq1, seq2):
    """Implements the Smith-Waterman local alignment algorithm."""
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    score_matrix = init_zero_matrix(rows, cols)

    # Fill in the matrix with local alignment scoring
    for i in range(1, rows):
        for j in range(1, cols):
            match = score_matrix[i-1][j-1] + score(i, j, seq1, seq2)
            delete = score_matrix[i-1][j] + GAP
            insert = score_matrix[i][j-1] + GAP
            score_matrix[i][j] = max(0, match, delete, insert)

    for row in score_matrix:
        print(row)  # optional: see matrix

    return traceback_sw(score_matrix, seq1, seq2)
