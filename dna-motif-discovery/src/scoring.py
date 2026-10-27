from collections import Counter

def extract_kmers(sequences, positions, k): # k-mers is a DNA seqence of length k, This function generates all possible k-mers 
    kmers = []
    for seq, pos in zip(sequences, positions): #loop through sequences and postions and pair them 
        kmers.append(seq[pos:pos + k]) #only append motif section to k-mers
    return kmers

def build_consensus(kmers): #Building consensus string by looking at each position and picking the most common letter
    
    k = len(kmers[0]) #length of motifs (all of them are the same length)
    consensus = ""
    
    for j in range(k): #go from 0 to length of motif
        column = [kmer[j] for kmer in kmers] # when j = 0, it will create a list of the first letters from every motif
        
        consensus += Counter(column).most_common(1)[0][0] # the counter functions counts how many times every letter has been present, then most_common letter is added to the consensus string
        
    return consensus
# THIS MAY NOT BE HAMMING DISTANCE DOUCLE CHECK 
def fitness_positions(sequences, positions, k):
    kmers = extract_kmers(sequences, positions, k) # extracts all k-mers
    consensus = build_consensus(kmers) #builds consensus motif from all k-mers
    score = 0
    for kmer in kmers: # goes through each kmer and compares it to consensus
        for a, b in zip(kmer, consensus): #basically hamming function 
            if a != b:
                score += 1 # if it is not the same it'll + 1 to fitness 

    return score, consensus, kmers
