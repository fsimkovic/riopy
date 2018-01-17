#!/usr/bin/env ccp4-python

from __future__ import print_function

import logging
import os
import sys
import time

from iotbx import pdb
from riopy.rio import Rio
from riopy.misc import calculate_origin_shift, move_model_to_target

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def select(h, sel):
    return h.select(h.atom_selection_cache().selection(sel))


def main(mtz, target, model, debug=False, dist=1.5, cells=1):

    if debug:
        logger.setLevel(logging.DEBUG)

    oshift = calculate_origin_shift(mtz, target, model)
    model = move_model_to_target(target, model, oshift)

    ref = pdb.pdb_input(target)
    uc = ref.crystal_symmetry().unit_cell()
    sg = ref.crystal_symmetry().space_group()

    h = ref.construct_hierarchy()
    s = pdb.pdb_input(model).construct_hierarchy()

    h_ca = select(h, "name ca")
    s_ca = select(s, "name ca")

    r = Rio(s_ca, h_ca, sg, uc)
    rscore = r.compute(max_dist=dist, cells=cells)
    logger.info("RIO score: %d", rscore)

    if debug:
        pass
    else:
        os.unlink(model)

    return rscore


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", default=False, action="store_true")
    parser.add_argument("-c", "--cells", default=1, type=int)
    parser.add_argument("-d", "--dist", default=1.5, type=float)
    parser.add_argument("mtz")
    parser.add_argument("target")
    parser.add_argument("model")
    main(**vars(parser.parse_args()))
