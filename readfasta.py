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

    seq1, seq2 = sequences[0], sequences[1] # Unpacking the sequences


    if is_dna(seq1) and is_dna(seq2) and user_input == 'dna':
        print("The sequences are DNA, as the user declared.")

    elif is_dna(seq1) and is_dna(seq2) and user_input == 'protein':
        print("The sequences are DNA, but the user declared them as protein. Check the input and try again.")

    elif is_protein(seq1) and is_protein(seq2) and user_input == 'protein':
        print("The sequences are protein, as the user declared.")
    elif is_protein(seq1) and is_protein(seq2) and user_input == 'dna':
        print("The sequences are protein, but the user declared them as DNA. Check the input and try again.")

    else:
        print("The sequences are either mixed or invalid. Check the input and try again")


# Ask the user for type of sequences: DNA or protein
user_input = get_sequence_type_from_user()
print(f"The user specified the sequences as: {user_input}")

# Reading the sequences from the FASTA file
seq1, seq2 = read_fasta("FASTA_diff_seq.txt")
print("Sequence 1:", seq1)
print("Sequence 2:", seq2)

sequences = [seq1, seq2] # Grouping the sequences into a list

# Checking if the sequences are of the same type as the user declared
sequence_type(sequences, user_input)