##################################################################################
######################### Smith-Waterman local alignment algorithm ###############
#################################################################################





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
        print("The sequences are DNA, but the user declared them as protein. Check the input and try again.")

    elif is_protein(sequences) and user_input == 'protein':
        print("The sequences are protein, as the user declared.")
    elif is_protein(sequences) and user_input == 'dna':
        print("The sequences are protein, but the user declared them as DNA. Check the input and try again.")

    else:
        print("The sequences are either mixed or invalid. Check the input and try again")


# Ask the user for type of sequences: DNA or protein
user_input = get_sequence_type_from_user()
print(f"The user specified the sequences as: {user_input}")

# Reading the sequences from the FASTA file
seq1, seq2 = read_fasta("FASTA_DNA.txt") # unpacking the sequences

print("Sequence 1:", seq1)
print("Sequence 2:", seq2)

sequences = [seq1, seq2] # Grouping the sequences into a list

# Checking if the sequences are of the same type as the user declared
sequence_type(sequences, user_input)

##################################################################################
############################## END OF FASTA HANDLING #############################
##################################################################################








##################################################################################
####################### 2. Smith-Waterman local alignment ###########################
##################################################################################

# Constants for scoring
MATCH = 1
MISMATCH = -1
GAP = -1


def score(i, j, seq1, seq2):
    """Simple scoring function: +1 for match, -1 for mismatch."""
    return MATCH if seq1[i-1] == seq2[j-1] else MISMATCH

def format_alignment(seq1, seq2):
    """Formats the aligned sequences above each other."""
    return f"{seq1}\n{seq2}"

# def init_zero_matrix(rows, cols):
#     """Initializes a score matrix filled with zeros for local alignment."""
#     return [[0 for _ in range(cols)] for _ in range(rows)]

# def reverse(seq1, seq2):
#     """Reverses two strings and returns them in a formatted string."""
#     #return f"Aligned Sequence 1: {seq1[::-1]}\nAligned Sequence 2: {seq2[::-1]}"

def traceback_sw(matrix, seq1, seq2):
    """Traceback for Smith-Waterman: starts at max value, ends at 0.
    Builds the aligned sequences from end to start."""
    # Finding max value in matrix
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

    while matrix[i][j] != 0:
        if i > 0 and j > 0 and matrix[i][j] == matrix[i-1][j-1] + score(i, j, seq1, seq2):
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif i > 0 and matrix[i][j] == matrix[i-1][j] + GAP:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1
        elif j > 0 and matrix[i][j] == matrix[i][j-1] + GAP:
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
        else:
            break
    return format_alignment(aligned_seq1, aligned_seq2)

    # Traceback until score is 0 to build the best local alignment
    # while matrix[i][j] != 0:
    #     if (i > 0 and j > 0 and 
    #         matrix[i][j] == matrix[i-1][j-1] + score(i, j, seq1, seq2)):
    #         aligned_seq1 += seq1[i-1] + aligned_seq1
    #         aligned_seq2 += seq2[j-1] + aligned_seq2
    #         i -= 1
    #         j -= 1
    #     elif i > 0 and matrix[i][j] == matrix[i-1][j] + GAP:
    #         aligned_seq1 += seq1[i-1] + aligned_seq1
    #         aligned_seq2 += "-" + aligned_seq2
    #         i -= 1
    #     else:
    #         aligned_seq1 += "-" + aligned_seq1
    #         aligned_seq2 += seq2[j-1] + aligned_seq2
    #         j -= 1

    # return reverse(aligned_seq2, aligned_seq1)  # reverse & return nicely



def smith_waterman(seq1, seq2):
    """Implements the Smith-Waterman local alignment algorithm and returns
    the best local alignment between two sequences."""
    rows = len(seq1) + 1
    cols = len(seq2) + 1
    #score_matrix = init_zero_matrix(rows, cols)
    score_matrix = [[0] * cols for _ in range(rows)]  # Initialize score matrix with zeros

    # Filling in the matrix with local alignment scoring
    for i in range(1, rows):
        for j in range(1, cols):
            match = score_matrix[i-1][j-1] + score(i, j, seq1, seq2)
            delete = score_matrix[i-1][j] + GAP
            insert = score_matrix[i][j-1] + GAP
            score_matrix[i][j] = max(0, match, delete, insert)

    return traceback_sw(score_matrix, seq1, seq2)

results = smith_waterman(seq1, seq2)
print(f"Smith-Waterman Local Alignment Result:\n{results}")