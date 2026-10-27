import random
from src.scoring import fitness_positions


def rmhc(sequences, k, iterations):
    """
    We represent a candidate solution as a list of start positions:
        positions[i] = start index of the k-mer chosen from sequences[i]

    Returns:
        best_positions, best_score, best_consensus
    """
    num_sequences = len(sequences)
    max_start_positions = [len(seq) - k for seq in sequences]

    # 1) Start with a random solution
    
    current_positions = [
        random.randint(0, max_start_positions[i])
        for i in range(num_sequences)
    ]
    current_score, current_consensus, _ = fitness_positions(
        sequences, current_positions, k
    )

    # Track the best solution seen so far (copy list with [:])
    best_positions = current_positions[:]
    best_score = current_score
    best_consensus = current_consensus

    # 2) Hill-climb for a number of iterations
    for _ in range(iterations):
        # Copy the current solution so we can "try" a mutation without
        # destroying the current solution unless it improves.
        new_positions = current_positions[:]

        # Pick ONE sequence index to mutate (small change in one dimension)
        i = random.randint(0, num_sequences - 1)

        # Small mutation: shift that one position by -1 or +1
        step = random.choice([-1, 1])
        new_pos = new_positions[i] + step

        # Clamp to valid range so the k-mer still fits inside the sequence
        # valid start positions are 0 .. (sequence_length - k)
        new_pos = max(0, min(max_start_positions[i], new_pos))
        new_positions[i] = new_pos

        # Score the mutated solution
        new_score, new_consensus, _ = fitness_positions(sequences, new_positions, k)

        # Greedy acceptance rule:
        # only accept the mutation if it IMPROVES (lowers) the score
        if new_score < current_score:
            current_positions = new_positions
            current_score = new_score
            current_consensus = new_consensus

            # Update best-ever (usually current == best in RMHC, but this is safe)
            if current_score < best_score:
                best_positions = current_positions[:]
                best_score = current_score
                best_consensus = current_consensus

    return best_positions, best_score, best_consensus