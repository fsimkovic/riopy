#!/usr/bin/env ccp4-python

import sys
import time

from iotbx import pdb
from riopy.rio import Rio


def select(h, sel):
    return h.select(h.atom_selection_cache().selection(sel))


def main():
    inputf = sys.argv[1]
    sminf = sys.argv[2]

    ref = pdb.pdb_input(inputf)
    uc = ref.crystal_symmetry().unit_cell()
    sg = ref.crystal_symmetry().space_group()

    h = ref.construct_hierarchy()
    s = pdb.pdb_input(sminf).construct_hierarchy()

    h_ca = select(h, "name ca")
    s_ca = select(s, "name ca")

    r = Rio(s_ca, h_ca, sg, uc)
    print(r.compute(max_dist=1.5, cells=1))


if __name__ == "__main__":
    main()
