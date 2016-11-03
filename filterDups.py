#!/usr/bin/env python2
"""Script for post processing duplication files by additional filtering based
on a given duplication size and percent identity
"""
from argparse import ArgumentParser


def filterDups(dupsf, dupSize, percentID):
    """Filters duplications from a dups file based on the supplied dupSize and
    percentID.
    """
    yield next(dupsf).rstrip()  # header
    for line in dupsf:
        stats = line.split()
        if int(stats[3]) >= dupSize and float(stats[4]) >= percentID:
            yield line.rstrip()


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('dups', help='dups file')
    parser.add_argument('-s', help='duplication size', type=int,
                        default=0, dest="dupSize")
    parser.add_argument('-%', help='percent identity', type=float,
                        default=0.0, dest="percentID")
    opts = parser.parse_args()

    with open(opts.dups, 'r') as dupsf:
        for dup in filterDups(dupsf, opts.dupSize, opts.percentID):
            print(dup)

if __name__ == '__main__':
    main()
