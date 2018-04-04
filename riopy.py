#!/usr/bin/env ccp4-python

from __future__ import print_function


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

    from riopy.rio import Rio
    r = Rio(args.mtz, args.target, args.model)
    score = r.compute(max_dist=args.dist, ncells=args.cells)
    print(score.__dict__)


if __name__ == "__main__":
    main()