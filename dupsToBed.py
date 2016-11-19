#!/usr/bin/env python2
"""Converts a given dups file from findScaffoldGaps into a bed file.
"""
import remove_overlap
from argparse import ArgumentParser


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('dups', help='dups file')
    opts = parser.parse_args()

    duplication = remove_overlap.parse_dups_file(opts.dups)
    for name, dups in duplication.items():
        for dup in dups:
            aboveStart = dup.gapStart - dup.dupSize
            behindEnd = dup.gapEnd + dup.dupSize
            aboveGap = "{}\t{}\t{}".format(name, aboveStart, dup.gapStart)
            belowGap = "{}\t{}\t{}".format(name, dup.gapEnd, behindEnd)
            print(aboveGap)
            print(belowGap)


if __name__ == '__main__':
    main()
