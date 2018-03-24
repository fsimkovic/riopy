import numpy as np
import os
import unittest

from riopy.rio import Rio


class TestRio(unittest.TestCase):
    def test_calculate_dists_2d_1(self):
        a = np.array([[1.000, 1.000, 1.000]])
        b = np.array([[1.000, 1.000, 1.000]])
        dists = Rio.calculate_dists_2d(a, b)
        self.assertListEqual([0.000], list(dists))

    def test_calculate_dists_2d_2(self):
        a = np.array([[1.000, 1.000, 1.000]])
        b = np.array([[-1.000, -1.000, -1.000]])
        dists = Rio.calculate_dists_2d(a, b)
        self.assertListEqual([3.464], list(dists))

    def test_calculate_dists_2d_3(self):
        a = np.array([[1.000, -1.000, 1.000]])
        b = np.array([[-1.000, 1.000, -1.000]])
        dists = Rio.calculate_dists_2d(a, b)
        self.assertListEqual([3.464], list(dists))

    def test_calculate_dists_2d_4(self):
        a = np.array([[1.000, -1.000, 1.000], [-1.000, 2.000, -1.000]])
        b = np.array([[-1.000, 1.000, -1.000], [2.000, -1.000, 2.000]])
        dists = Rio.calculate_dists_2d(a, b)
        self.assertListEqual([3.464, 5.196], list(dists))

    def test_calculate_dists_2d_5(self):
        a = np.array([[4.543, -16.459, 19.186], [3.859, -20.237, 18.659], [6.162, -20.346, 15.576]])
        b = np.array([[3.916, -16.906, 18.009], [3.708, -20.627, 18.760], [5.258, -21.443, 15.386]])
        dists = Rio.calculate_dists_2d(a, b)
        self.assertListEqual([1.407, 0.430, 1.434], list(dists))

    def test_find_contacts_1(self):
        self.assertListEqual([], list(Rio.find_contacts(np.array([]), 0.0)))

    def test_find_contacts_2(self):
        self.assertListEqual([], list(Rio.find_contacts(np.array([1.000]), 0.0)))

    def test_find_contacts_3(self):
        self.assertListEqual([1], list(Rio.find_contacts(np.array([1.000, 0.500, 0.510]), 0.5)))

    def test_find_contacts_4(self):
        self.assertListEqual([0, 1], list(Rio.find_contacts(np.array([-1.000, 0.500, 0.510]), 0.5)))


if __name__ == "__main__":
    unittest.main(verbosity=2)
