import random
from src.scoring import fitness_positions


def random_search(sequences, k, iterations):
    """
    Random search over position guesses.

    Returns:
        best_positions, best_score, best_consensus
    """
    num_sequences = len(sequences)

    # Maximum valid start index for each sequence
    max_start_positions = [len(seq) - k for seq in sequences]

    best_score = float("inf")  # set the best score to infinity
    best_positions = None
    best_consensus = None

    for _ in range(iterations):
        # 1. random guess: one position per sequence
        positions = [
            random.randint(0, max_start_positions[i])
            for i in range(num_sequences)
        ]

        # 2. score it
        score, consensus, _ = fitness_positions(sequences, positions, k)

        # 3 keep best score so far
        if score < best_score:
            best_score = score
            best_positions = positions
            best_consensus = consensus

    return best_positions, best_score, best_consensus