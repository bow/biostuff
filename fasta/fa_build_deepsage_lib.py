#!/usr/bin/env python

"""
Build SAGE tag library from Ensembl genes for mapping with Bowtie 1.0.0.


This script creates a FASTA file containing all SAGE tags from all Ensembl
transcript (ENST) records. Identical tags from the same gene (occurs when a
gene has multiple isoforms with the same 3'-most tag) are collapsed into a
single record. The script does not disambiguate identical tags from different
genes.

Required files:

    * An ENST-ENSG mapping file
      Two column tab-delimited file that maps an ENST ID to an ENSG ID

    * FASTA file containing ENST records from UCSC
      The script assumes that all FASTA IDs are formatted according to UCSC's
      FASTA formatting. If your FASTA file comes from another source, you will
      need to change the way the script extracts ENST ID from each FASTA record.

Requirements:
    * Python == 2.7.x
    * Biopython >= 1.6.0 <http://biopython.org>

Copyright (c) 2013 Wibowo Arindrarto <w.arindrarto@lumc.nl>
Copyright (c) 2013 LUMC Sequencing Analysis Support Core <sasc@lumc.nl>
MIT License <http://opensource.org/licenses/MIT>

"""

RELEASE = False
__version_info__ = ('0', '1', )
__version__ = '.'.join(__version_info__)
__version__ += '-dev' if not RELEASE else ''


import argparse
import re

from Bio import SeqIO


# regex pattern for extracting ENST IDs
_ID_PAT = re.compile(r'hg19_.*_(ENST\d+)')
# regex pattern for extracting SAGE tags; assumes no ambiguous bases
_TAG_PAT = re.compile(r'CATG[CATG]{17,17}?')


def parse_map(map_file):
    """Given a path to an ENST-ENSG map file, return a dictionary of all
    mappings.

    :param map_file: path to map file
    :type map_file: str
    :returns: ENST-ENSG ID mapping
    :rtype: dict

    """
    with open(map_file, 'r') as src:
        lines = (line.strip().split() for line in src)
        return {key: value for key, value in lines}


def find_tags(rec):
    """Given a SeqRecord object, return the 3'-most SAGE tag.

    :param rec: input sequence record
    :type rec: SeqRecord
    :returns: 3'-most SAGE tag
    :rtype: str

    """
    matches = re.findall(_TAG_PAT, str(rec.seq))
    if matches:
        return matches[-1]


def build_tag_lib(transcriptome_file, map_file):
    """Builds a SAGE tag library for each ENSG record.

    :param transcriptome_file: path to FASTA file containing ENST records
    :type transcriptome_file: str
    :param map_file: path to ENST-ENSG mapping file
    :type map_file: str
    :returns: a dictionary of tags for each ENSG record, keys are ENSG ID and
              values are lists of tag sequences for the record.
    :rtype: dict {str: list}

    """
    tag_lib = {}
    ens_map = parse_map(map_file)
    for rec in SeqIO.parse(transcriptome_file, 'fasta'):
        enst_id = re.search(_ID_PAT, rec.id).group(1)
        ensg_id = ens_map[enst_id]
        tag = find_tags(rec)
        if tag:
            if ensg_id not in tag_lib:
                # using set, so identical tags are counted once
                tag_lib[ensg_id] = set([tag])
            else:
                tag_lib[ensg_id].add(tag)

    return tag_lib


if __name__ == '__main__':

    usage = __doc__.split('\n\n\n')
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=usage[0], epilog=usage[1])

    parser.add_argument('transcriptome', type=str,
            help='Path to FASTA transcriptome file')
    parser.add_argument('id_map', type=str,
            help='Path to ENST-ENSG map file')
    parser.add_argument('tags', type=str,
            help='Path to output FASTA tag file')
    parser.add_argument('--version', action='version', version='%(prog)s ' +
            __version__)

    args = parser.parse_args()

    tag_lib = build_tag_lib(args.transcriptome, args.id_map)

    with open(args.tags, 'w') as fasta:
        # sort tag library based on how many tags an ENSG record has
        # (descending)
        tag_sorted = sorted(tag_lib.items(), key=lambda x: len(x[1]), 
                reverse=True)
        for gene, tags in tag_sorted:
            for idx, tag in enumerate(tags):
                # for each ENSG tag, append '_{num}' to the ENSG ID
                tag_id = '%s_%s' % (gene, str(idx + 1))
                fasta.write(">%s\n%s\n" % (tag_id, tag))
