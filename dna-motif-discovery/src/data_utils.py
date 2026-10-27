import random

DNA = "ACGT"


# Here we are defining a random sequence of DNA of n length 
def rand_dna(n):
    result = ""
    for i in range(n):
        letter = random.choice(DNA) # the .choice function is a varaible from random module, it picks a item at random from the sequence
        result += letter
    return result

def hamming(a, b):
    count = 0
    for x, y in zip(a, b):  # pairs each item of the same index of different strings together, this line unpacks into x & y varaibles 
        if x != y:
            count += 1
    return count



def make_noisy_copy(motif, max_mismatches):
    motif_list = list(motif) ## Change it into a list so we can change the motif  e.g "ACTG' -> ['A', 'C', 'T', 'G']
    
    # mismatches must be shorter than the actual length of the motif 
    max_mismatches = min(max_mismatches, len(motif))
    
    # Generate random no. of mismatches from 0 to max_mismatches (which is determined by the user)
    k = random.randint(0, max_mismatches)
    
    # choose which positions to change
    positions_to_change = random.sample(range(len(motif)), k) # e.g. if motif length = 6 and k = 2, we could get [0,3] or [4,5]
    
    for idx in positions_to_change:  # the index here refers to the positions that will be changed for. e.g. [0,3] first idx = 0 and second idx = 3 
        original = motif_list[idx] # here we call for the original letter 
        
        # append, for loop which has a if loop  
        # here we append each letter (a/g/c/t) to the 'choices', only if the letter is not the same one 
        choices = [letter for letter in DNA if letter != original] 
      
        motif_list[idx] = random.choice(choices)
        
    return "".join(motif_list)


def generate_dataset(motif, num_sequences, sequence_length, max_mismatches):
    sequences = []
    true_positions = []
    planted_copies = []

    for i in range(num_sequences):
        seq = rand_dna(sequence_length)
        pos = random.randint(0, sequence_length - len(motif))

        planted = make_noisy_copy(motif, max_mismatches)
        new_seq = seq[:pos] + planted + seq[pos + len(motif):]

        sequences.append(new_seq)
        true_positions.append(pos)
        planted_copies.append(planted)

    return sequences, true_positions, planted_copies