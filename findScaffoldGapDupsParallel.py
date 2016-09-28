#!/usr/bin/env python
"""Find the weird tandem duplications around scaffold gaps and report
their location and size.

Requires sonLib, jobTree, and the pypi package twobitreader."""
import re
import os
from argparse import ArgumentParser
from collections import namedtuple
from twobitreader import TwoBitFile
from sonLib.bioio import fastaRead, popenCatch, getTempFile
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack

# The HoxD70 scoring matrix, modified so that matches with Ns have
# constant -120 score.
hoxd70 = { 'A': { 'A':   91, 'C': -114, 'G':  -31, 'T': -123, 'N': -120 },
           'C': { 'A': -114, 'C':  100, 'G': -125, 'T':  -31, 'N': -120 },
           'G': { 'A':  -31, 'C': -125, 'G':  100, 'T': -114, 'N': -120 },
           'T': { 'A': -123, 'C':  -31, 'G': -114, 'T':   91, 'N': -120 },
           'N': { 'A': -120, 'C': -120, 'G': -120, 'T': -120, 'N': -120 },
   }

def alignWithConstrainedLeftEdge(seq1, seq2, xDrop=200):
    """Gapless extension from left edge with x-drop filtering.

    Returns size, score,  and % ID of the match.

    This code would be useful if the dups *both* started at the
    start/end of the scaffold gap, but since one starts and one ends
    at the scaffold gap, this is totally useless. Keeping it as a
    monument to my stupidity.
    """
    peakSize = 0
    peakScore = 0
    peakID = 0
    curSize = 0
    curScore = 0
    curID = 0

    # Trim sequences to be the same length
    if len(seq1) > len(seq2):
        seq1 = seq1[:len(seq2)]
    if len(seq2) > len(seq1):
        seq2 = seq2[:len(seq1)]

    # Do gapless extension
    for c1, c2 in zip(seq1, seq2):
        curSize += 1
        curScore += hoxd70[c1][c2]
        curID += 1 if c1 == c2 else 0
        if peakScore - curScore >= xDrop:
            break
        elif curScore >= peakScore:
            peakScore = curScore
            peakSize = curSize
            peakID = curID

    if peakSize > 0:
        percentID = float(peakID)/peakSize
    else:
        percentID = 0.0
    return peakSize, peakScore, percentID

def alignWithLastz(seq1, seq2):
    """Special-case of aligning two sequences for which the first part of
    seq2 may be similar to the last part of seq1.

    Returns size and % ID of the alignment."""
    seq2Path = getTempFile()
    try:
        with open(seq2Path, 'w') as seq2File:
            seq2File.write('>\n' + seq2)
        mappingStr = popenCatch("lastz --strand=plus --step=4 --maxwordcount=10 %s[unmask]"
                                " --nogapped --xdrop=300 --identity=95 "
                                " --format=mapping-" % seq2Path, stdinString=">\n" + seq1.upper())
        # There may be multiple alignments, but there can be at most one
        # that fits our requirements (starts at 0 in seq1, ends at
        # len(seq2))

        for line in mappingStr.split("\n"):
            if len(line) == 0:
                continue
            fields = line.split()
            # idPct looks looks like "100%": remove the % and convert to
            # float
            identity = float(fields[8][:-1])
            # Starts are 0-based inclusive, ends are 0-based exclusive
            seq2Start = int(fields[1])
            seq2End = int(fields[2])
            seq1Start = int(fields[5])
            seq1End = int(fields[6])
            if seq1End != len(seq1) or seq2Start != 0:
                # Invalid.
                continue
            # If we get here, this is the one possible alignment, so
            # it's safe to return without checking the rest
            return seq2End - seq2Start, identity
        return 0, 0.0
    finally:
        os.remove(seq2Path)

def findGaps(sequence):
    for match in re.finditer(r'[Nn]+', sequence):
        yield (match.start(), match.end())

Gap = namedtuple('Gap', ['header', 'start', 'end', 'before', 'after'])

class Concatenate(Target):
    def __init__(self, inputs, output):
        Target.__init__(self)
        self.inputs = inputs
        self.output = output

    def run(self):
        with open(self.output, 'w') as outfile:
            # Print header
            outfile.write('sequence\tgapStart\tgapEnd\tdupSize\tdupPctID\n')
            for input in self.inputs:
                with open(input) as infile:
                    for line in infile:
                        outfile.write(line)

class AlignAndCompare(Target):
    def __init__(self, twoBit, gaps, output):
        Target.__init__(self, memory=2000000000)
        self.twoBit = twoBit
        self.gaps = gaps
        self.output = output

    def run(self):
        genome = TwoBitFile(self.twoBit)
        with open(self.output, 'w') as outfile:
            for gap in self.gaps:
                print gap
                print repr(genome[gap.header][gap.before:gap.start].upper())
                print repr(genome[gap.header][gap.end:gap.after].upper())
                size, percentID = alignWithLastz(genome[gap.header][gap.before:gap.start].upper(), genome[gap.header][gap.end:gap.after].upper())
                if size > 20:
                    outfile.write('\t'.join(map(str, [gap.header, gap.start, gap.end, size, percentID])) + '\n')

class FindScaffoldGapsForTwoBit(Target):
    def __init__(self, twoBit, output, maxSize, split):
        Target.__init__(self, memory=4000000000)
        self.twoBit = twoBit
        self.output = output
        self.maxSize = maxSize
        self.split = split

    def run(self):
        gaps = []
        for header, sequence in TwoBitFile(self.twoBit).iteritems():
            for gapStart, gapEnd in findGaps(str(sequence)):
                beforeGap = max(0, gapStart - self.maxSize)
                afterGap = min(len(sequence), gapEnd + self.maxSize)
                if gapStart - beforeGap > 5 and afterGap - gapEnd > 5:
                    gaps.append(Gap(header=header, start=gapStart, end=gapEnd, before=beforeGap, after=afterGap))

        outputs = []
        for gapList in [gaps[i:i + self.maxSize] for i in xrange(0, len(gaps), self.maxSize)]:
            output = getTempFile(rootDir=self.getGlobalTempDir())
            self.addChildTarget(AlignAndCompare(self.twoBit, gapList, output))
            outputs.append(output)

        self.setFollowOnTarget(Concatenate(outputs, self.output))

class FindScaffoldGapsForAllTwoBits(Target):
    def __init__(self, twoBits, outputs, maxSize, split):
        Target.__init__(self)
        self.twoBits = twoBits
        self.outputs = outputs
        self.maxSize = maxSize
        self.split = split

    def run(self):
        for twoBit, output in zip(self.twoBits, self.outputs):
            self.addChildTarget(FindScaffoldGapsForTwoBit(twoBit, output, self.maxSize, self.split))

def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--twoBits', help='2bit files', required=True, nargs='+')
    parser.add_argument('--outputs', help='output files', required=True, nargs='+')
    parser.add_argument('--maxSize', help='maximum size to attempt to check', type=int,
                        default=5000)
    parser.add_argument('--split', help='Allow this many checks per job', type=int,
                        default=1000)
    Stack.addJobTreeOptions(parser)
    opts = parser.parse_args()

    Stack(FindScaffoldGapsForAllTwoBits(opts.twoBits, opts.outputs, opts.maxSize, opts.split)).startJobTree(opts)

if __name__ == '__main__':
    from findScaffoldGapDupsParallel import *
    main()
