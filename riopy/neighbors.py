
__author__ = "Felix Simkovic"
__date__ = "16 Jan 2018"
__version__ = "1.0"

import itertools

from halocell import HaloCell
from symops import SymmetryOperator
from scitbx.array_family import flex


class Neighbors(object):

    @staticmethod
    def ops(sg, n_uc):
        """Yield all rotation / translation operations in cartesian space
        
        Parameters
        ----------
        sg : str
           The space group
        n_uc : int
           The number of neighboring cells

        Yields
        ------
        list
           A generator of (rotation, translation) operators

        """
        for symop, cellop in itertools.product(SymmetryOperator.ops(sg), HaloCell.ops(n_uc)):
            rotop = symop.r().as_double()
            transop = tuple(flex.double(symop.t().as_double()) + cellop)
            yield (rotop, transop)

