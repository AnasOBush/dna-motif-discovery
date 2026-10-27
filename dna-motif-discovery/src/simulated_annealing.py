import random
import math
from src.scoring import fitness_positions


def simulated_annealing(
    sequences,
    k,
    iterations=5000,
    initial_temperature=10.0,
    cooling_rate=0.995,
    min_temperature=0.001
):
    """
    Simulated Annealing for DNA motif discovery.

    A solution is represented as a list of motif start positions:
        positions[i] = start index of the k-mer chosen from sequences[i]

    The algorithm:
    - starts from a random solution
    - proposes a small mutation
    - always accepts better moves
    - sometimes accepts worse moves with a probability depending on temperature
    - gradually lowers the temperature over time

    Returns:
        best_positions, best_score, best_consensus
    """

    if any(len(seq) < k for seq in sequences):
        raise ValueError("All sequences must have length at least k.")

    num_sequences = len(sequences)

    # Maximum valid start index for each sequence
    max_start_positions = [len(seq) - k for seq in sequences]

    # 1) Start with a random solution
    current_positions = [
        random.randint(0, max_start_positions[i])
        for i in range(num_sequences)
    ]

    current_score, current_consensus, _ = fitness_positions(
        sequences, current_positions, k
    )

    # Track best solution found
    best_positions = current_positions[:]
    best_score = current_score
    best_consensus = current_consensus

    temperature = initial_temperature

    # 2) Main loop
    for _ in range(iterations):
        # Stop early if temperature becomes too small
        if temperature < min_temperature:
            break

        # Copy current solution
        new_positions = current_positions[:]

        # Mutate one randomly chosen sequence position
        i = random.randint(0, num_sequences - 1)

        # Small mutation: shift by -1 or +1
        step = random.choice([-1, 1])
        new_pos = new_positions[i] + step

        # Clamp to valid range for this specific sequence
        new_pos = max(0, min(max_start_positions[i], new_pos))
        new_positions[i] = new_pos

        # Evaluate new solution
        new_score, new_consensus, _ = fitness_positions(sequences, new_positions, k)

        # Difference in score (lower is better)
        delta = new_score - current_score

        # Accept if better, or sometimes if worse
        if delta < 0:
            accept = True
        else:
            probability = math.exp(-delta / temperature)
            accept = random.random() < probability

        if accept:
            current_positions = new_positions
            current_score = new_score
            current_consensus = new_consensus

            if current_score < best_score:
                best_positions = current_positions[:]
                best_score = current_score
                best_consensus = current_consensus

        # Cool down
        temperature *= cooling_rate

    return best_positions, best_score, best_consensus