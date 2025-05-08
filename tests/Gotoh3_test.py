import numpy as np


##################################################################################
###################### 1. READING AND HANDLING FASTA FILES ##########################
##################################################################################

def read_fasta(file):
    '''Reads a FASTA file and returns the sequences as strings in a list.'''
    sequences = [] # empty list to store sequences
    with open(file, 'r') as file:
        seq = "" # empty string
        for line in file: # read line by line
            if line.startswith('>'): # if it's a header line
                if seq: # if there's a sequence already
                    sequences.append(seq) # append to list
                seq = ""
            else: # if not a header line
                seq += line.strip().upper() # remove whitespace and convert to uppercase
        if seq: # if there's a sequence left at the end
            sequences.append(seq) # append the last sequence

    if len(sequences) != 2:
        raise ValueError("FASTA file must contain exactly two sequences.")
    return sequences[0], sequences[1] # return the two sequences

def is_dna(sequences):
    '''Checks if the sequences are DNA.'''
    valid_nucleotides = {'A', 'C', 'G', 'T', 'N'}
    for seq in sequences:
        if not all(nucleotide in valid_nucleotides for nucleotide in seq):
            return False
    return True


def is_protein(sequences):
    '''Checks if the sequences are protein'''
    valid_aa = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}
    dna_only = {'A', 'C', 'G', 'T'}
    for seq in sequences:
        # Checks if all characters are valid amino acids, otherwise returns False
        if not all(aa in valid_aa for aa in seq):
            return False
        
        # Checks if all characters are A, G, C, T, then it's likely DNA
        if all(aa in dna_only for aa in seq):
            return False
    return True


def get_sequence_type_from_user():
    '''Asks the user to specify if the sequences are DNA or protein.'''
    while True:
        user_input = input("Are the sequences DNA or protein? (Enter 'DNA' or 'protein'):").strip().lower()
        if user_input in {'dna', 'protein'}:
            return user_input
        else: 
            print("Invalid input. Please enter 'DNA' or 'protein'.")


def sequence_type(sequences, user_input):
    '''Checks if the two sequences are of the same type, as the user says (either DNA or protein).'''

    if is_dna(sequences) and user_input == 'dna':
        print("The sequences are DNA, as the user declared.")

    elif is_dna(sequences) and user_input == 'protein':
        raise TypeError("The sequences are DNA, but the user declared them as protein. Check the input and try again.")
    
    elif is_protein(sequences) and user_input == 'protein':
        print("The sequences are protein, as the user declared.")
    elif is_protein(sequences) and user_input == 'dna':
        raise TypeError("The sequences are protein, but the user declared them as DNA. Check the input and try again.")

    else:
        raise TypeError("The sequences are either mixed or invalid. Check the input and try again")

# Ask the user for type of sequences: DNA or protein
user_input = get_sequence_type_from_user()
print(f"The user specified the sequences as: {user_input}")

# Reading the sequences from the FASTA file
seq1, seq2 = read_fasta('Pairwise-alignment-project\\tests\\test_fasta') # unpacking the sequences

print("Sequence 1:", seq1)
print("Sequence 2:", seq2)

sequences = [seq1, seq2] # Grouping the sequences into a list

# Checking if the sequences are of the same type as the user declared
sequence_type(sequences, user_input)


##################################################################################
############################## END OF FASTA HANDLING #############################
##################################################################################




##################################################################################
####################### 2. Needleman-Wunsch global alignment ###########################
##################################################################################

# Constants for scoring
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

# A, B = read_fasta('Pairwise-alignment-project\\tests\\test_fasta')
alignment_with_prints(seq1,seq2)
# alignment_with_prints("CG", "CCGA")
# alignment_with_prints("GATTACA", "CATGCATCGAC")
