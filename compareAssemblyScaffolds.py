#!/usr/bin/env python2
"""For a given assembly in twobit file format, as well its corresponding
duplications file, returns the assemblies total number of scaffold gaps,
the number of scaffold gaps found to be apart of a duplication, and a
percentage based on the ratio of these two values.
"""
from argparse import ArgumentParser
from twobitreader import TwoBitFile
from findScaffoldGapDups import findGaps


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('twobit', help='twobit assembly file')
    parser.add_argument('dups', help='corresponding dups file')
    opts = parser.parse_args()

    with open(opts.dups, 'r') as dupsf:
        next(dupsf)  # Header
        dupsScaffolds = sum(1 for line in dupsf)

    total = 0
    for sequence in TwoBitFile(opts.twobit).itervalues():
        total += sum(1 for _, _ in findGaps(str(sequence)))

    base = opts.dups.split('/')[-1]
    base = base.split('.')[0]
    print("{}\t{}\t{}\t{:.2f}".format(base, total, dupsScaffolds,
          (float(dupsScaffolds) / total) * 100))


if __name__ == '__main__':
    main()
