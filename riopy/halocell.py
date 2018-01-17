
import itertools

from scitbx.array_family import flex


class HaloCell(object):

    @staticmethod
    def ops(ndim):
        mates = set(((0, 0, 0), ))
        def add_one_layer():
            operations = [0, 1, -1]
            for mate in mates.copy():
                for x_op, y_op, z_op in itertools.product([0, 1, -1], repeat=3):
                    new_mate = (mate[0] + x_op, mate[1] + y_op, mate[2] + z_op)
                    mates.add(new_mate)
        for _ in range(ndim):
            add_one_layer()
        return [flex.double(m) for m in mates]
