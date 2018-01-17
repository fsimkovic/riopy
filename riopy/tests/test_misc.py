import unittest

from riopy.misc import *


class Test(unittest.TestCase):
    def test_count_res_in_frag_1(self):
        nfrag = count_res_in_frag([1, 2, 3], min_frag_len=3)
        self.assertEqual(3, nfrag)

    def test_count_res_in_frag_2(self):
        nfrag = count_res_in_frag([], min_frag_len=3)
        self.assertEqual(0, nfrag)

    def test_count_res_in_frag_3(self):
        nfrag = count_res_in_frag([1, 2, 3, 11, 12, 13], min_frag_len=3)
        self.assertEqual(6, nfrag)

    def test_count_res_in_frag_4(self):
        nfrag = count_res_in_frag([1, 2, 3, 11, 12, 13, 14, 15, 21, 22], min_frag_len=3)
        self.assertEqual(8, nfrag)

    def test_count_res_in_frag_5(self):
        nfrag = count_res_in_frag([1, 2, 3, 11, 12, 13, 14, 15, 21, 22], min_frag_len=5)
        self.assertEqual(5, nfrag)

    def test_count_res_in_frag_6(self):
        nfrag = count_res_in_frag([1], min_frag_len=2)
        self.assertEqual(0, nfrag)

    def test_count_res_in_frag_7(self):
        res = [83, 87, 104, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122]
        nfrag = count_res_in_frag(res, min_frag_len=3)
        self.assertEqual(16, nfrag)


if __name__ == "__main__":
    unittest.main(verbosity=2)
