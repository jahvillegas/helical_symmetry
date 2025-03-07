import unittest
from src.core.helical_symmetry import compute_helical_parameters

class TestHelicalSymmetry(unittest.TestCase):

    def test_compute_helical_parameters(self):
        matrix, vector = compute_helical_parameters(10, 5, 30, 4.8)
        self.assertEqual(vector.tolist(), [0, 0, 4.8])
        self.assertAlmostEqual(matrix[0][0], 0.866, places=3)

if __name__ == "__main__":
    unittest.main()

