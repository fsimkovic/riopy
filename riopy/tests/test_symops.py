
import unittest

from riopy.symops import SymmetryOperators

class SymmetryOperatorsTest(unittest.TestCase):

    def test___init___1(self):
        symops = SymmetryOperators.all_ops("P1")
        self.assertTrue(len(symops) == 1)
        self.assertTupleEqual((0.0, 0.0, 0.0), symops[0].t().as_double())
        self.assertTupleEqual((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0), symops[0].r().as_double())




if __name__ == "__main__":
    unittest.main(verbosity=2)
