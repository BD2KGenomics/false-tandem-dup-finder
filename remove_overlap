#!/usr/bin/env python
"""
Remove overlap (+ extra if desired) from tandem duplications
around scaffold gaps.
"""
import sys
from argparse import ArgumentParser
from collections import defaultdict, namedtuple

from sonLib.bioio import fastaRead, fastaWrite

FalseDup = namedtuple('FalseDup', ['gapStart', 'gapEnd', 'dupSize'])

def parse_dups_file(path):
    """
    Parse a .dups file produced by findScaffoldGapDups.
    """
    dups = defaultdict(list)
    with open(path) as f:
        header = f.readline().strip().split()
        if header != ['sequence', 'gapStart', 'gapEnd', 'dupSize', 'dupPctID']:
            raise RuntimeError("Unexpected .dups file header: %s" % header)
        for line in f:
            fields = line.strip().split()
            dups[fields[0]].append(FalseDup(int(fields[1]), int(fields[2]), int(fields[3])))
    return dups

def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('fasta', help='Sequence file')
    parser.add_argument('dups', help='File describing the scaffold gap dups')
    parser.add_argument('--additional', type=int, default=0,
                        help='Additional amount to trim past the dup')
    opts = parser.parse_args()
    # Ingest dup locations
    dups = parse_dups_file(opts.dups)
    for header, sequence in fastaRead(opts.fasta):
        sequence_dups = dups[header]
        offset = 0
        sequence = list(sequence)
        for dup in sequence_dups:
            del sequence[dup.gapEnd - offset:dup.gapEnd + dup.dupSize + opts.additional - offset]
            del sequence[dup.gapStart - dup.dupSize - opts.additional - offset:dup.gapStart - offset]
            offset += 2 * (dup.dupSize + opts.additional)
        fastaWrite(sys.stdout, header, ''.join(sequence))

if __name__ == '__main__':
    main()
