#!/usr/bin/env python2
"""Find the weird tandem duplications around scaffold gaps and report
their location and size.

Requires sonLib, jobTree, and the pypi package twobitreader."""
from argparse import ArgumentParser
from collections import namedtuple
from twobitreader import TwoBitFile
from sonLib.bioio import getTempFile
from jobTree.scriptTree.target import Target
from jobTree.scriptTree.stack import Stack
from findScaffoldGapDups import alignWithBlat, findGaps


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
                print(gap)
                print(repr(genome[gap.header][gap.before:gap.start].upper()))
                print(repr(genome[gap.header][gap.end:gap.after].upper()))
                seq1 = genome[gap.header][gap.before:gap.start].upper()
                seq2 = genome[gap.header][gap.end:gap.after].upper()
                size, percentID = alignWithBlat(seq1, seq2)
                if size > 20:
                    stats = [gap.header, gap.start, gap.end, size, percentID]
                    outfile.write(('\t'.join(str(s) for s in stats)) + '\n')


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
                    gaps.append(Gap(header=header, start=gapStart, end=gapEnd,
                                    before=beforeGap, after=afterGap))

        outputs = []
        for i in xrange(0, len(gaps), self.maxSize):
            for gapList in [gaps[i:i + self.maxSize]]:
                output = getTempFile(rootDir=self.getGlobalTempDir())
                align = AlignAndCompare(self.twoBit, gapList, output)
                self.addChildTarget(align)
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
            find = FindScaffoldGapsForTwoBit(twoBit, output, self.maxSize,
                                             self.split)
            self.addChildTarget(find)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--twoBits', help='2bit files', required=True,
                        nargs='+')
    parser.add_argument('--outputs', help='output files', required=True,
                        nargs='+')
    parser.add_argument('--maxSize', help='maximum size to attempt to check',
                        type=int, default=5000)
    parser.add_argument('--split', help='Allow this many checks per job',
                        type=int, default=1000)
    Stack.addJobTreeOptions(parser)
    opts = parser.parse_args()

    finds = FindScaffoldGapsForAllTwoBits(opts.twoBits, opts.outputs,
                                          opts.maxSize, opts.split)
    Stack(finds).startJobTree(opts)


if __name__ == '__main__':
    from findScaffoldGapDupsParallel import *
    main()
