import random
from src.scoring import fitness_positions


def _random_chromosome(max_start_positions):
    """
    A chromosome = one position per sequence.
    Example (5 sequences): [3, 10, 7, 0, 21]
    """
    return [
        random.randint(0, max_start_positions[i])
        for i in range(len(max_start_positions))
    ]


def _tournament_selection(population, scores, tournament_size=3):
    """
    Pick 'tournament_size' random individuals and return the best (lowest score).
    """
    n = len(population)
    best_idx = None

    for _ in range(tournament_size):
        idx = random.randint(0, n - 1)
        if best_idx is None or scores[idx] < scores[best_idx]:
            best_idx = idx

    return population[best_idx]


def _one_point_crossover(parent_a, parent_b):
    """
    One-point crossover:
    - take first part from parent_a
    - take second part from parent_b
    """
    n = len(parent_a)
    if n < 2:
        return parent_a[:], parent_b[:]

    cut = random.randint(1, n - 1)
    child1 = parent_a[:cut] + parent_b[cut:]
    child2 = parent_b[:cut] + parent_a[cut:]
    return child1, child2


def _mutate(chromosome, max_start_positions, mutation_rate=0.1):
    """
    Mutate each gene with probability 'mutation_rate'.
    Small mutation: shift by -1 or +1 (like RMHC).
    """
    child = chromosome[:]

    for i in range(len(child)):
        if random.random() < mutation_rate:
            step = random.choice([-1, 1])
            new_pos = child[i] + step
            new_pos = max(0, min(max_start_positions[i], new_pos))
            child[i] = new_pos

    return child


def genetic_algorithm(sequences, k, population_size=30, generations=200, mutation_rate=0.1, tournament_size=3):
    """
    Genetic Algorithm over position lists.

    Returns:
        best_positions, best_score, best_consensus
    """
    # Maximum valid start index for each sequence
    max_start_positions = [len(seq) - k for seq in sequences]

    # 1) initial population (random chromosomes)
    population = [
        _random_chromosome(max_start_positions)
        for _ in range(population_size)
    ]

    # Track best ever
    best_positions = None
    best_score = float("inf")
    best_consensus = None

    for _ in range(generations):
        # 2) evaluate population
        scores = []
        consensuses = []
        for chrom in population:
            s, cons, _ = fitness_positions(sequences, chrom, k)
            scores.append(s)
            consensuses.append(cons)

            if s < best_score:
                best_score = s
                best_positions = chrom[:]
                best_consensus = cons

        # 3) build next generation (elitism: keep best 1)
        # find current generation best
        gen_best_idx = min(range(len(scores)), key=lambda i: scores[i])
        elite = population[gen_best_idx][:]

        new_population = [elite]

        # 4) fill the rest using selection + crossover + mutation
        while len(new_population) < population_size:
            parent_a = _tournament_selection(population, scores, tournament_size)
            parent_b = _tournament_selection(population, scores, tournament_size)

            child1, child2 = _one_point_crossover(parent_a, parent_b)

            child1 = _mutate(child1, max_start_positions, mutation_rate)
            child2 = _mutate(child2, max_start_positions, mutation_rate)

            new_population.append(child1)
            if len(new_population) < population_size:
                new_population.append(child2)

        population = new_population

    return best_positions, best_score, best_consensus