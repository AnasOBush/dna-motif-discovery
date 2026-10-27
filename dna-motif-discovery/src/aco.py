import random
from src.scoring import fitness_positions


def aco_motif_search(
    sequences,
    k,
    num_ants=20,
    iterations=100,
    evaporation_rate=0.2,
    alpha=1.0,
    q=1.0
):
    """
    Ant Colony Optimisation for motif discovery using position representation.

    Returns:
        best_positions, best_score, best_consensus
    """
    num_sequences = len(sequences)

    # Maximum valid start index for each sequence
    max_start_positions = [len(seq) - k for seq in sequences]

    # Pheromone[seq_idx][start_pos]
    # Each sequence may have a different number of valid start positions
    pheromone = [
        [1.0 for _ in range(max_start_positions[seq_idx] + 1)]
        for seq_idx in range(num_sequences)
    ]

    best_positions = None
    best_score = float("inf")
    best_consensus = None

    for _ in range(iterations):
        ant_solutions = []

        # Each ant builds one candidate solution
        for _ in range(num_ants):
            positions = []

            for seq_idx in range(num_sequences):
                num_positions = max_start_positions[seq_idx] + 1

                weights = [
                    pheromone[seq_idx][pos] ** alpha
                    for pos in range(num_positions)
                ]

                chosen_pos = random.choices(
                    range(num_positions),
                    weights=weights,
                    k=1
                )[0]

                positions.append(chosen_pos)

            score, consensus, _ = fitness_positions(sequences, positions, k)
            ant_solutions.append((positions, score, consensus))

            if score < best_score:
                best_positions = positions[:]
                best_score = score
                best_consensus = consensus

        # Evaporation
        for seq_idx in range(num_sequences):
            num_positions = max_start_positions[seq_idx] + 1

            for pos in range(num_positions):
                pheromone[seq_idx][pos] *= (1.0 - evaporation_rate)

                # avoid pheromone becoming too tiny
                if pheromone[seq_idx][pos] < 0.0001:
                    pheromone[seq_idx][pos] = 0.0001

        # Deposit pheromone from all ants, stronger for better solutions
        for positions, score, _ in ant_solutions:
            deposit = q / (1.0 + score)

            for seq_idx, pos in enumerate(positions):
                pheromone[seq_idx][pos] += deposit

    return best_positions, best_score, best_consensus