#!/usr/bin/env python2
"""Find the weird tandem duplications around scaffold gaps and report
their location and size."""
import re
import os
# import subprocess
from argparse import ArgumentParser
from sonLib.bioio import fastaRead, system, getTempFile


def alignWithBlat(seq1, seq2):
    """Special-case of aligning two sequences for which the first part of
    seq2 may be similar to the last part of seq1.

    Returns size and % ID of the alignment."""
    seq1Path = getTempFile()
    seq2Path = getTempFile()
    outPath = getTempFile()
    try:
        with open(seq1Path, 'w') as seq1File:
            with open(seq2Path, 'w') as seq2File:
                seq1File.write('>\n' + seq1)
                seq2File.write('>\n' + seq2)

        # Run blat command
        # system("blat %s %s -q=dna -tileSize=12 -stepSize=4 -minIdentity=95"
        #            " -repMatch=10 -noHead %s" % (seq1Path, seq2Path, outPath))
        cmd = "blat {} {} -q=dna -tileSize=12 -stepSize=4 -minIdentity=95" \
              " -repMatch=10 -noHead -fastMap {}".format(seq1Path, seq2Path, outPath)
        # with open(os.devnull, 'w') as devnull:
        #     subprocess.call(cmd, shell=True, bufsize=-1, stdout=devnull)
        system(cmd)


        # There may be multiple alignments, but there can be at most one
        # that fits our requirements (starts at 0 in seq1, ends at
        # len(seq2))

        with open(outPath, 'r') as outFile:
            for line in outFile:
                line = line.rstrip()
                # print(line)
                if len(line) == 0:
                    continue
                stats = line.split('\t')
                # Starts are 0-based inclusive, ends are 0-based exclusive
                seq2Start = int(stats[11])
                seq1End = int(stats[16])
                if seq1End != len(seq1) or seq2Start != 0:
                    # Invalid.
                    continue
                # If we get here, this is the one possible alignment, so
                # it's safe to return without checking the rest
                # return seq2End - seq2Start, identity

                dupPct = float(stats[0]) / (float(stats[0]) + float(stats[1]))
                dupPct *= 100
                # TODO num of mismatches or size of duplication?
                return stats[0], dupPct
            return 0, 0.0
    finally:
        os.remove(seq1Path)
        os.remove(seq2Path)
        os.remove(outPath)


def findGaps(sequence):
    """Generator yielding the start and ends of a scaffold with any number
    of Ns.
    """
    patt = re.compile(r"[Nn]+")
    for match in patt.finditer(sequence):
        yield (match.start(), match.end())


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('fasta', help='fasta file')
    parser.add_argument('--maxSize', help='maximum size to attempt to check',
                        default=5000)
    opts = parser.parse_args()

    # Print header
    print('sequence\tgapStart\tgapEnd\tdupSize\tdupPctID')

    for header, seq in fastaRead(opts.fasta):
        for gapStart, gapEnd in findGaps(seq):
            beforeGap = seq[max(0, gapStart - opts.maxSize):gapStart].upper()
            afterGap = seq[gapEnd:min(len(seq), gapEnd + opts.maxSize)].upper()
            if len(beforeGap) < 5 or len(afterGap) < 5:
                # N at the end/beginning of a sequence, without enough
                # information on one or both ends to be useful
                continue
            size, percentID = alignWithBlat(beforeGap, afterGap)
            if size > 20:
                read_info = [header, gapStart, gapEnd, size, percentID]
                format_str = "{}\t{}\t{}\t{}\t{:.1f}"
                print(format_str.format(*read_info))

if __name__ == '__main__':
    main()
