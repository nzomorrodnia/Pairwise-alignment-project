#!/home/tgc/anaconda3/bin/python3

import numpy as np
import textwrap as tw
import sys 

##################################################################################
###################### 1. READING AND HANDLING FASTA FILES ##########################
##################################################################################

def read_fasta(file):
    '''Reads a FASTA file and returns the sequences as strings in a list.'''
    try:
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
    except IOError as error:
        return f"Could not open file because: {error}"



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
seq1, seq2 = read_fasta(fr'{sys.argv[1]}') # unpacking the sequences


print("Sequence 1:", seq1)
print("Sequence 2:", seq2)

sequences = [seq1, seq2] # Grouping the sequences into a list

# Checking if the sequences are of the same type as the user declared
sequence_type(sequences, user_input)

##################################################################################
############################## END OF FASTA HANDLING #############################
##################################################################################




##################################################################################
####################### 2. Needleman-Wunsch global alignment #####################
##################################################################################

# Constants for scoring
MATCH = 1
MISMATCH = -1
GAP_OPEN = -10
GAP_EXTEND = -1

# Setting custom score parameters
custom_constants = input('Do you want custom score parameters? Type Y if "Yes":\n').lower()
if custom_constants == 'y':
    terminate = False
    while not terminate:
        try:
            MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND = (map(int, user_constants:=input('Please input MATCH, MISMATCH, GAP_OPEN and GAP_EXTEND as integers separated by spaces, otherwise type ENTER for default parameters:\n').split()))
            if GAP_EXTEND < GAP_OPEN:
                print('Gap extension cannot be more negative than gap opening')
                continue
            if MATCH < 0:
                print('Match must be non-negative')
                continue
            if MISMATCH > 0:
                print('Mismatch must be zero or negative')
                continue
            terminate = True
        except ValueError:
            if not bool(user_constants):
                terminate = True
                break
            else:
                print('Please enter exactly four integers')
                print(user_constants)
                continue


def get_matrices(n, m):
    '''Initializes matrices'''
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
    '''Backtracking through the three matrices. The function will decide which matrix to use 
    based on the available paths for the current matrix and call itself until 'i' and 'j' are 0,
    and then append current path to a list'''
    if i == 0 and j == 0:
        alignments.append(current_alignment)
        return

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
    '''Takes the aligned strings and reverses them'''
    alignments_rev = []
    for alignment in alignments:
        l_seq1 = list(alignment[0])
        l_seq1.reverse()
        l_seq2 = list(alignment[1])
        l_seq2.reverse()
        alignments_rev.append((''.join(l_seq1), ''.join(l_seq2)))

    return list(dict.fromkeys(alignments_rev))


def needleman_wunsch(seq1, seq2):
    '''Main alignment function. Fills initialized matrices from get_matrices(), calls
    traceback() and sends output to all_reversed_and_no_dupes()'''
    n, m = len(seq1), len(seq2)
    S, P, Q = get_matrices(n, m)

    # Filling initialized matrices based on applied scoring rules
    for i in range(1, n+1):
        for j in range(1, m+1):
            match_score = MATCH if seq1[i - 1] == seq2[j - 1] else MISMATCH
            P[i][j] = max(S[i-1][j] + GAP_OPEN + GAP_EXTEND, P[i - 1][j] + GAP_EXTEND)
            Q[i][j] = max(S[i][j-1] + GAP_OPEN + GAP_EXTEND, Q[i][j - 1] + GAP_EXTEND)
            S[i][j] = max(S[i-1][j-1] + match_score, P[i][j], Q[i][j])

    # Traceback starting point
    max_score = max(S[n][m], P[n][m], Q[n][m])
    alignments = []
    if S[n][m] == max_score:
        traceback(seq1, seq2, S, P, Q, n, m, ('', ''), alignments, 'S')
    if P[n][m] == max_score:
        traceback(seq1, seq2, S, P, Q, n, m, ('', ''), alignments, 'P')
    if Q[n][m] == max_score:
        traceback(seq1, seq2, S, P, Q, n, m, ('', ''), alignments, 'Q')

    return all_reversed_and_no_dupes(alignments)


def alignment_with_prints(seq1, seq2):
    '''Prints alignments from previous outputs'''
    raw_alignments = needleman_wunsch(seq1, seq2)
    
    # The function output will be processed to be properly displayed on screen
    result_list = list()
    sep = 40    # Separator for screen output

    for res in raw_alignments:
        # Collecting symbols to display
        symbols = ""
        for n in range(len(res[0])):
            if res[0][n] == res[1][n]:
                symbols += "|"
            elif res[0][n] == '-':
                symbols += ' '
            elif res[1][n] == '-':
                symbols += ' '
            else:
                symbols += '*' 
        
        # Wrapping the two aligned sequences and its symbols to fit the screen
        wp_res0 = tw.wrap(res[0], width=sep)
        wp_res1 = tw.wrap(res[1], width=sep)
        symbols = tw.wrap(symbols, width=sep)

        # Collecting sequence positions for each sequence    
        seq1_posmarks = [(pos := pos + sep - wp_res0[n - 1].count('-') if n > 0 else 1) for n in range(len(res[0]) // sep + 1)]
        seq2_posmarks = [(pos := pos + sep - wp_res1[n - 1].count('-') if n > 0 else 1) for n in range(len(res[1]) // sep + 1)]

        # Collecting all data in a dictionary for systematic accessing
        result = {
            'seq1_ALIGNED': {'slices':wp_res0,'slice_lengths':seq1_posmarks},
            'seq2_ALIGNED': {'slices':wp_res1,'slice_lengths':seq2_posmarks},
            'sym_list': symbols
        }
        result_list.append(result)


    # Final print statement
    print(f"for:\n\n{ seq1}\n\nand:\n\n{seq2}\n\ngot:\n")
    for result in result_list: 
        fragments_to_display = result['seq1_ALIGNED']['slices']
        for n in range(len(fragments_to_display)):
            print(f"{result['seq1_ALIGNED']['slice_lengths'][n]}\t"
                  f"{result['seq1_ALIGNED']['slices'][n]}"
                  f"\n\t{result['sym_list'][n]}"
                  f"\n{result['seq2_ALIGNED']['slice_lengths'][n]}"
                  f"\t{result['seq2_ALIGNED']['slices'][n]}\n")
        print('---------------------------------')
    
    print(f'Finished, {len(raw_alignments)} possible alignments found.Scores used:\n'
          f'Match = {MATCH}, Mismatch = {MISMATCH}, Gap opening = {GAP_OPEN}, Gap extension = {GAP_EXTEND}')



alignment_with_prints(seq1,seq2)
