#!/usr/bin/env python

"""
Quick and dirty script for removing N's from FASTQ records.


By default only Ns at the 3' ends are removed. To remove Ns
in the 5', use the --remove-5prime-n flag.

Requirements:
    * Python == 2.7.x
    * Biopython >= 1.63

Author: w.arindrarto@lumc.nl
"""

import argparse
import re

from Bio import SeqIO

_RE_END_N = re.compile(r'[Nn]+$')
_RE_START_N = re.compile(r'^[Nn]+')


def trim_ns(infile, outfile, encoding, min_length, remove_5prime):

    def trim(infile):
        for rec in SeqIO.parse(infile, 'fastq-' + encoding):
            ends = _RE_END_N.search(str(rec.seq))
            if ends:
                rec = rec[:-len(_RE_END_N.search(str(rec.seq)).group(0))]

            if remove_5prime:
                starts = _RE_START_N.search(str(rec.seq))
                if starts:
                    rec = rec[len(_RE_END_N.search(str(rec.seq)).group(0)):]

            if len(rec) >= min_length:
                yield rec

    trimmed = trim(infile)

    with open(outfile, 'w') as dest:
        written = SeqIO.write(trimmed, dest, 'fastq-' + encoding)

    print "Written", written, "FASTQ records."


if __name__ == '__main__':

    usage = __doc__.split('\n\n\n')
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=usage[0], epilog=usage[1])

    parser.add_argument('input', type=str, help='Input FASTQ file')
    parser.add_argument('output', type=str, help='Output FASTQ file')
    parser.add_argument('--encoding', type=str, choices=['sanger', 'illumina', 'solexa'],
        default='sanger', help='FASTQ encoding [default: sanger]')   
    parser.add_argument('--min-length', type=int, default=20,
        help='Minimum length of FASTQ record to keep')   
    parser.add_argument('--remove-5prime-n', dest='remove_5prime_n', default=False,
        action='store_true', help='Whether to remove 5\' Ns or not [default: False]')

    args = parser.parse_args()

    trim_ns(args.input, args.output, args.encoding, args.min_length, args.remove_5prime_n)
