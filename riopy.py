#!/usr/bin/env ccp4-python

from __future__ import print_function

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", default=False, action="store_true")
    parser.add_argument("-c", "--cells", default=1, type=int)
    parser.add_argument("-d", "--dist", default=1.5, type=float)
    parser.add_argument("mtz")
    parser.add_argument("target")
    parser.add_argument("model")
    args = parser.parse_args()

    purge = True
    if args.debug:
        logger.setLevel(logging.DEBUG)
        purge = False 

    from riopy.rio import Rio
    r = Rio(args.mtz, args.target, args.model, purge=purge)
    score = r.compute(max_dist=args.dist, ncells=args.cells)
    logger.info("Overall RIO: %.3f", score.total)
    logger.info("RIO model normalised: %.3f", score.norm_model)
    logger.info("RIO target normalised: %.3f", score.norm_target)


if __name__ == "__main__":
    main()
