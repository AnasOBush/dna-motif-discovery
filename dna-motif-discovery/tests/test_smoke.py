"""Small deterministic sanity checks for the custom motif-search implementations."""

import random
import unittest

from src.data_utils import generate_dataset, hamming
from src.scoring import fitness_positions
from src.random_search import random_search
from src.rmhc import rmhc
from src.ga import genetic_algorithm
from src.aco import aco_motif_search
from src.simulated_annealing import simulated_annealing


class MotifTests(unittest.TestCase):
    def test_hamming(self):
        self.assertEqual(hamming("ACGT", "ACGA"), 1)
        self.assertEqual(hamming("ACGT", "ACGT"), 0)

    def test_synthetic_generator(self):
        random.seed(17)
        seqs, positions, planted = generate_dataset("ACG", 5, 12, 0)
        self.assertEqual(len(seqs), 5)
        for seq, pos, motif in zip(seqs, positions, planted):
            self.assertEqual(seq[pos:pos + 3], motif)
            self.assertEqual(motif, "ACG")

    def test_fitness(self):
        score, consensus, kmers = fitness_positions(["TTACGG", "GGACGT"], [2, 2], 3)
        self.assertEqual(score, 0)
        self.assertEqual(consensus, "ACG")
        self.assertEqual(kmers, ["ACG", "ACG"])

    def test_searches_return_valid_candidates(self):
        random.seed(17)
        sequences = ["GGACGTT", "TTACGTA", "CCACGTC"]
        algorithms = (
            lambda: random_search(sequences, 3, iterations=20),
            lambda: rmhc(sequences, 3, iterations=20),
            lambda: genetic_algorithm(sequences, 3, population_size=8, generations=5),
            lambda: aco_motif_search(sequences, 3, num_ants=5, iterations=5),
            lambda: simulated_annealing(sequences, 3, iterations=20),
        )
        for algo in algorithms:
            positions, score, consensus = algo()
            self.assertEqual(len(positions), len(sequences))
            self.assertEqual(len(consensus), 3)
            self.assertGreaterEqual(score, 0)
            self.assertTrue(all(0 <= pos <= len(seq) - 3 for pos, seq in zip(positions, sequences)))


if __name__ == "__main__":
    unittest.main()
