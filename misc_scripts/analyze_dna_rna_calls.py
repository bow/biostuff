#!/usr/bin/env python2.7

"""
Allelic expression finder.


This script compares calls from at least two VCF files: one from DNA-seq and
the other from RNA-seq to find candidate positions where expression is
monoallelic.

Monoallelicly-expressed regions are deduced from positions where the RNA-seq
variant is homozygous (reference or alternative) and the DNA-seq variant is
heterozygous.

At the moment, this script is designed to work with RNA-seq variant calls from
VarScan and DNA-seq variant calls from samtools. Each variant file must only
have one sample calls. RNA-seq VCF file must be created using VarScan's
mpileup2cns, which lists variant as well as reference calls. DNA-seq VCF file
must be annotated using vcf-annotate (from the vcftools suite), otherwise
the gene columns will all be '?'.

Requirements:
    * Python == 2.7.x
    * PyVCF > 0.6.4

Copyright (c) 2013 Wibowo Arindrarto <w.arindrarto@lumc.nl>
Copyright (c) 2013 LUMC Sequencing Analysis Support Core <sasc@lumc.nl>
MIT License <http://opensource.org/licenses/MIT>
"""

from __future__ import print_function

RELEASE = False
__version_info__ = ('0', '1', )
__version__ = '.'.join(__version_info__)
__version__ += '-dev' if not RELEASE else ''


import argparse
import os
import sys
from collections import namedtuple, OrderedDict

import vcf
from vcf import utils

# column order of per sample results files
PERSAMPLE_COLS = [
    'CHROM', 'POS', 'GENE', 'REF_GT', 'RNA_ZYGOS',
    'DNA_GT', 'DNA_GQ', 'DNA_DEPTH', 'DNA_DEPTH_REF/ALT',
    'RNA_GT', 'RNA_GQ', 'RNA_DEPTH', 'RNA_DEPTH_REF/ALT',
    'FOUND_IN'
]
# line format of per sample result files
PERSAMPLE_LINEFMT = '{' + '}\t{'.join(PERSAMPLE_COLS) + '}'
# column order of overlap file
OVERLAP_COLS = [
    'CHROM', 'POS', 'GENE', 'REF_GT', 'RNA_ZYGOS',
    'DNA_GT', 'DNA_GQ', 'DNA_DEPTH', 'DNA_DEPTH_REF/ALT',
    'RNA_GTS', 'RNA_GQS', 'RNA_DEPTHS', 'RNA_DEPTHS_REF/ALT'
]
# line format of overlap file
OVERLAP_LINEFMT = '{' + '}\t{'.join(OVERLAP_COLS) + '}'


# simple container for DNA-seq calls
DnaCall = namedtuple('DnaCall', [
    'chrom', 'pos', 'gt_ref', 'gt', 'gq', 'depth', 'depth_ref', 'depth_alt',
])

def _build_dna_call(call):
    """Returns a DnaCall namedtuple given a VCF record object from a samtools VCF ."""
    gt_ref = call.samples[0].gt_phase_char().join(str(call.alleles[x]) for
            x in (0, 0))
    gt = call.samples[0].gt_bases
    gq = call.samples[0].data.GQ
    depth = sum([(int(x)) for x in call.INFO['DP4']])
    depth_ref = sum([(int(x)) for x in call.INFO['DP4'][:2]])
    depth_alt = sum([(int(x)) for x in call.INFO['DP4'][2:4]])

    return DnaCall(call.CHROM, call.POS, gt_ref, gt, gq, depth, depth_ref,
            depth_alt)


def _vars_type(rna_calls):
    """Given a list of VCF record objects from a VarScan output,
    returns a string denoting the zygosity composition of the calls."""
    het_count, hom_count = 0, 0
    for rna_call in rna_calls:
        if rna_call.samples[0].is_het:
            het_count += 1
        else:
            hom_count += 1

    if het_count == 0:
        assert hom_count != 0
        return 'ALL_HOM'
    elif hom_count == 0:
        assert het_count != 0
        return 'ALL_HET'
    return 'MIX'


class Results(object):

    "Class representing DNA & RNA variant intersection results."""

    def __init__(self, out_dir, dna_vcf, *rna_vcfs):

        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        self.out_dir = out_dir

        self.dna_reader = vcf.Reader(filename=dna_vcf)
        self.dna_name = os.path.basename(dna_vcf)

        rna_readers = []
        rna_names = []
        outs_per_sample = OrderedDict()
        for rna_vcf in rna_vcfs:
            basename = os.path.basename(rna_vcf)
            rna_names.append(basename)
            rna_readers.append(vcf.Reader(filename=rna_vcf))
            # also create open file handles for each sample to write results to
            out_file = os.path.join(out_dir, os.path.splitext(basename)[0])
            out_file = out_file + '.results'
            outf = open(out_file, 'w')
            outf.write('\t'.join(PERSAMPLE_COLS) + '\n')
            outs_per_sample[basename] = outf

        self.rna_readers = rna_readers
        self.rna_names = rna_names
        self.outs_per_sample= outs_per_sample

        outov = open(os.path.join(out_dir, 'overlap.results'), 'w')
        outov.write('#RNA_SAMPLE_ORDER: ' + ', '.join(rna_names) + '\n')
        outov.write('\t'.join(OVERLAP_COLS) + '\n')
        self.outs_overlap = outov

        samerec = lambda rec: (rec.CHROM, rec.POS, rec.REF)

        for calls in utils.walk_together(self.dna_reader, *self.rna_readers,
                                         vcf_record_sort_key=samerec):
            # select for heterozygous calls present in DNA and any one of the RNAs
            if calls[0] is not None and calls[0].samples[0].is_het and \
                    any(calls[1:]):
                # gene annotation
                gene = '?'
                if 'ANN' in calls[0].INFO:
                    gene = ','.join(calls[0].INFO['ANN'])

                if all(calls[1:]):
                    self.write_overlap(gene, calls[0], calls[1:])
                    print('Found in DNA and all RNA:', calls[0].CHROM,
                          calls[0].POS, gene, file=sys.stderr)
                else:
                    print('Found in DNA and {0} RNA:'.format(len([x for x in
                          calls[1:] if x])), calls[0].CHROM, calls[0].POS,
                          gene, file=sys.stderr)

                self.write_per_rna(gene, calls[0], calls[1:])

    def write_overlap(self, gene, dna_call, *rna_calls):
        dna_call = _build_dna_call(dna_call)
        cols = {
            'CHROM': dna_call.chrom,
            'POS': dna_call.pos,
            'GENE': gene,
            'REF_GT': dna_call.gt_ref,
            'DNA_GT': dna_call.gt,
            'DNA_GQ': dna_call.gq,
            'DNA_DEPTH': dna_call.depth,
            'DNA_DEPTH_REF/ALT': '{0}/{1}'.format(dna_call.depth_ref,
                dna_call.depth_alt),
            'RNA_GTS': '',
            'RNA_GQS': '',
            'RNA_DEPTHS': '',
            'RNA_DEPTHS_REF/ALT': '',
        }

        # iterate over each RNA call, processing ones that are not None
        for idx, rna_call in enumerate(*rna_calls):
            refd = int(rna_call.samples[0].data.RD)
            altd = int(rna_call.samples[0].data.AD)
            cols['RNA_GTS'] += ',' + str(rna_call.samples[0].gt_bases)
            cols['RNA_GQS'] += ',' + str(rna_call.samples[0].data.GQ)
            cols['RNA_DEPTHS'] += ',' + str(refd + altd)
            cols['RNA_DEPTHS_REF/ALT'] += ',' + '{0}/{1}'.format(refd, altd)
            cols['RNA_ZYGOS'] = _vars_type(*rna_calls)
        # remove start commas
        cols['RNA_GTS'] = cols['RNA_GTS'][1:]
        cols['RNA_GQS'] = cols['RNA_GQS'][1:]
        cols['RNA_DEPTHS'] = cols['RNA_DEPTHS'][1:]
        cols['RNA_DEPTHS_REF/ALT'] = cols['RNA_DEPTHS_REF/ALT'][1:]

        self.outs_overlap.write(OVERLAP_LINEFMT.format(**cols) + '\n')

    def write_per_rna(self, gene, dna_call, *rna_calls):
        dna_call = _build_dna_call(dna_call)

        found_in = []
        for idx, rna_call in enumerate(*rna_calls):
            if rna_call is not None and not rna_call.samples[0].is_het:
                found_in.append(self.rna_names[idx])
        found_in = ','.join(found_in)

        # iterate over each RNA call, processing ones that are not None
        for idx, rna_call in enumerate(*rna_calls):
            # skipping None calls
            if rna_call is None:
                continue

            refd = int(rna_call.samples[0].data.RD)
            altd = int(rna_call.samples[0].data.AD)
            cols = {
                'CHROM': dna_call.chrom,
                'POS': dna_call.pos,
                'GENE': gene,
                'REF_GT': dna_call.gt_ref,
                'DNA_GT': dna_call.gt,
                'DNA_GQ': dna_call.gq,
                'DNA_DEPTH': dna_call.depth,
                'DNA_DEPTH_REF/ALT': '{0}/{1}'.format(dna_call.depth_ref,
                    dna_call.depth_alt),
                'RNA_GT': rna_call.samples[0].gt_bases,
                'RNA_GQ': rna_call.samples[0].data.GQ,
                'RNA_DEPTH': refd + altd,
                'RNA_DEPTH_REF/ALT': '{0}/{1}'.format(refd, altd),
                'FOUND_IN': found_in,
            }
            if rna_call.samples[0].is_het:
                cols['RNA_ZYGOS'] = 'HET'
            else:
                cols['RNA_ZYGOS'] = 'HOM'
            # the list order are maintained, that's why we can do this
            # to fetch which file object to write to
            rna_out = self.outs_per_sample[cur_rna_name]
            rna_out.write(PERSAMPLE_LINEFMT.format(**cols) + '\n')


if __name__ == '__main__':

    usage = __doc__.split('\n\n\n')
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=usage[0], epilog=usage[1])

    parser.add_argument('--dna', type=str,
            help='Path to VCF file of DNA variant calls')
    parser.add_argument('--rna', nargs='+', type=str,
            help='Path to VCF file of RNA-seq variant calls')
    parser.add_argument('--out-dir', type=str, default='allelic_expression',
            help='Path to results output directory')
    parser.add_argument('--version', action='version', version='%(prog)s ' +
            __version__)

    args = parser.parse_args()

    Results(args.out_dir, args.dna, *args.rna)
