#!/usr/bin/env python2

"""
Spliting a single BAM file into multiple BAM files based on the reference.


For a given BAM file, this script will produce multiple BAM files, each
containing reads mapped to a single record in the BAM reference. For example, a
BAM file from a full human genome sequencing will be split into at least 24-25
files (autosomes, sex chromosomes, mitochondrial genome).

By default, the script writes all split BAMs into the current directory.
For unmapped reads, you may choose whether to output them into a separate file
or not.

The splitting may be sped up by indexing the BAM file beforehand using samtools
index. The index must be located in the same directory as the input BAM file,
with the '.bam.bai' extension. Note that this will not allow you to output
unmapped BAM reads into its own file, however.

Requirements:
    * Python == 2.7.x
    * pysam >= 0.7 <http://www.cgat.org/~andreas/documentation/pysam/api.html>

Copyright (c) 2013 Wibowo Arindrarto <w.arindrarto@lumc.nl>
Copyright (c) 2013 LUMC Sequencing Analysis Support Core <sasc@lumc.nl>
MIT License <http://opensource.org/licenses/MIT>

"""

RELEASE = False
__version_info__ = ('0', '1', )
__version__ = '.'.join(__version_info__)
__version__ += '-dev' if not RELEASE else ''


import argparse
import itertools
import os
import warnings

import pysam


def split_by_fetch(input_bam, out_dir):
    """Splits a BAM file into multiple BAM files each with reads mapped to the
    same reference sequence by fetching the reads.

    Requires that the BAM file is indexed.

    :param input_bam: path to input BAM file
    :type input_bam: str
    :param out_dir: path to output directory containing all split BAMs
    :type out_dir: str
    :returns: None

    """
    input_bam_index = input_bam + '.bai'
    assert os.path.exists(input_bam_index), "BAM index file %r not found." \
            % input_bam_index

    in_bam = pysam.Samfile(input_bam, 'rb')

    for ref in in_bam.references:
        ref_records = in_bam.fetch(ref)

        # peek into first record, so we don't have to write header-only BAM
        # files later on
        try:
            first_rec = ref_records.next()
        except StopIteration:
            continue

        out_path = os.path.join(out_dir, ref + '.bam')
        out_bam = pysam.Samfile(out_path, 'wb', template=in_bam)
        print "Writing to %r" % out_path
        for rec in itertools.chain([first_rec], ref_records):
            out_bam.write(rec)
        out_bam.close()


def split_by_iter(input_bam, out_dir, unmapped_bam=None):
    """Splits a BAM file into multiple BAM files each with reads mapped to the
    same reference sequence by iterating over all reads.

    :param input_bam: path to input BAM file
    :type input_bam: str
    :param out_dir: path to output directory containing all split BAMs
    :type out_dir: str
    :param unmapped_bam: file name of BAM file for unmapped reads
    :type unmapped_bam: str
    :returns: None

    """
    in_bam = pysam.Samfile(input_bam, 'rb')
    if unmapped_bam is not None:
        un_path = os.path.join(out_dir, unmapped_bam)
        if os.path.exists(un_path):
            raise ValueError("File %r already present" % un_path)
        un_bam = pysam.Samfile(un_path, 'wb', template=in_bam)
    else:
        un_bam = None
    cur_ref, cur_bam = None, None

    un_count = 0
    for rec in in_bam:
        if rec.tid < 0:
            # any records with tid == 0 should be unmapped
            assert rec.flag & 0x4
            un_count += 1
            if un_bam is not None:
                un_bam.write(rec)
            continue
        # get actual record name
        rec_ref = in_bam.getrname(rec.tid)

        # if we encounter a new tid, open new file and write to it
        if cur_ref != rec_ref or (cur_ref is None and cur_bam is None):
            cur_ref = rec_ref
            cur_path = os.path.join(out_dir, cur_ref + '.bam')
            # TODO: allow for non-reference-sorted BAM files? For now, just
            # raise an error in for those files
            if os.path.exists(cur_path):
                raise ValueError("File %r already present" % cur_path)
            if cur_bam is not None:
                cur_bam.close()
            print "Writing to %r" % cur_path
            cur_bam = pysam.Samfile(cur_path, 'wb', template=in_bam)

        cur_bam.write(rec)

    print "Found %i unmapped records" % un_count
    if un_bam is not None:
        print "Unmapped records have been written to %r" % un_path


if __name__ == '__main__':

    usage = __doc__.split('\n\n\n')
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=usage[0], epilog=usage[1])

    parser.add_argument('input_bam', type=str,
            help='Input BAM file')
    parser.add_argument('-d', '--output-directory', dest='out_dir',
            type=str, default=os.getcwd(),
            help='Output directory')
    parser.add_argument('-i', '--index', dest='index',
            default=False, action='store_true',
            help='Whether to use BAM index or not')
    parser.add_argument('-u', '--unmapped-bam', dest='unmapped_bam',
            default=None, help='Path to write unmapped BAM records')
    parser.add_argument('--version', action='version', version='%(prog)s ' +
            __version__)

    args = parser.parse_args()

    # check for input file names and output directory presence
    assert os.path.exists(args.input_bam), "Input BAM file %r does not exist" % \
            args.input_bam
    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    if args.index is None:
        split_by_iter(args.input_bam, args.out_dir, args.unmapped_bam)
    else:
        if args.unmapped_bam is not None:
            warnings.warn("Supplied unmapped BAM file name %r will not be used for"
                " splitting using BAM indexes." % args.unmapped_bam)
        split_by_fetch(args.input_bam, args.out_dir)
