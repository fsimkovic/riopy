
from cctbx.crystal import symmetry

class SymmetryOperator(object):

    @staticmethod
    def ops(sg):
        return symmetry(space_group=sg).space_group().all_ops()

