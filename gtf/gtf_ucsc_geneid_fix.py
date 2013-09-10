#!/usr/bin/env python2.7

"""
Updates 'gene_id' entries in a GTF file downloaded from UCSC Table Browser
to have gene IDs as values instead of transcript IDs.


Two types annotation sources can be used to replace the 'gene_id' values:

1. Local annotation source. To use this, supply a file name to the
   '--local' argument. The file must only have two columns denoting
   the transcript - gene ID mapping (the first column contain the
   transcript IDs).

2. Remote annotation source (UCSC). To use this, supply the UCSC
   database to use (e.g. 'hg19') to the '--db' argument and the annotation
   source to the '--annot' argument. Annotation source is either 'ucsc',
   'ensembl', 'refseq', 'gencode14', or 'gencode17'.

You can only choose local or remote sources, not both.

Requirements:
    * Python == 2.7.x
    * MySQLdb >= 1.2.3
    * track >= 1.3.0-dev (dev version from github.com/xapple/track)

Copyright (c) 2013 Wibowo Arindrarto <w.arindrarto@lumc.nl>
Copyright (c) 2013 LUMC Sequencing Analysis Support Core <sasc@lumc.nl>
MIT License <http://opensource.org/licenses/MIT>

"""

RELEASE = False
__version_info__ = ('0', '1', )
__version__ = '.'.join(__version_info__)
__version__ += '-dev' if not RELEASE else ''


import argparse
import os
import warnings

import track
import MySQLdb


# Credentials for connecting to database
CRED = {
    'host': 'genome-mysql.cse.ucsc.edu',
    'user': 'genome',
}

# Queries that return ('transcript ID', 'gene ID') on various tables
QUERIES = {
    'ucsc': 'SELECT knownGene.name, kgXref.geneSymbol FROM ' \
            'knownGene, kgXref WHERE knownGene.name = kgID',
    'ensembl': 'SELECT name, name2 FROM ensGene',
    'refseq': 'SELECT name, name2 FROM refGene',
    'gencode17': 'SELECT name, name2 FROM wgEncodeGencodeBasicV17',
    'gencode14': 'SELECT name, name2 FROM wgEncodeGencodeBasicV14',
}


def get_ucsc_transcript_gene_mapping(annot_src, db, cred=CRED):
    """Returns the transcript-gene name mapping for an annotation source
    from a given database source.

    The function does not know whether the annotation source exists within
    the given database nor does it do any check before trying connect.

    :param annot_src: name of annotation source
    :type annot_src: str
    :param db: UCSC database name to use
    :param cred: database login credentials, must contain entries for at
            least 'host' and 'user', defaults to credentials for
            public UCSC server
    :type cred: dict
    :returns: transcript-gene name mapping
    :rtype: dict(str: str)

    """
    con = MySQLdb.connect(db=db, **cred)
    cur = con.cursor()
    cur.execute(QUERIES[annot_src])

    return {tid: gid for tid, gid in cur.fetchall()}


def get_local_transcript_gene_mapping(fname):
    """Loads a two-column file (transcript ID - gene ID) as a dict.

    :param fname: path to file
    :type fname: str
    :returns: mapping of column 1 and column 2 in the file
    :rtype: dict(str: str)

    """
    mapping = {}
    with open(fname, 'r') as src:
        for line in src:
            if not line.strip():
                break
            elif not line:
                continue

            key, value = line.strip().split()
            if key in mapping:
                if value != mapping[key]:
                    raise ValueError("Duplicate transcript ID ({0}) with "
                            "ambiguous gene ID ({1} vs {2}).".format(
                            key, value, mapping[key]))
            mapping[key] = value

    return mapping


def update_gene_id_attr(chrom_recs, mapping):
    """Given an iterator for `track` records, yield `track` records with
    updated gene ID.

    :param chrom_recs: iterator returning `track` records for one chromosome
    :type chrom_recs: iterator
    :returns: generator yielding single `track` records
    :rtype: (yield) `track.pyrow.SuperRow`

    """
    for rec in chrom_recs:
        data = list(rec.data)
        # gene ID is always at index 7 (first index of attributes)
        init_gid = data[7]
        try:
            map_gid = mapping[init_gid]
        except KeyError:
            warnings.warn("Transcript ID {0} not found in the given "
                "mapping, initial value is left unchanged.".format(
                init_gid))
        else:
            data[7] = map_gid

        yield data


def ucsc_geneid_fix(in_gtf, out_gtf, remote=None, local=None):
    """Updates 'gene_id' entries in GTF files downloaded from UCSC
    Table Browser to contain gene IDs instead of transcript IDs.

    If the output GTF file name already exists, it will be overwritten.

    :param in_gtf: path to input GTF file
    :type in_gtf: str
    :param out_gtf: path to output GTF file
    :type out_gtf: str
    :param remote: UCSC database and annotation source to use
    :type remote: dict('db': str, 'annot_src': str)
    :param local: two-column file name containing transcript-gene mapping,
            only when `db` and `annot_src` are None
    :type local: str
    :returns: None

    """
    # remote not defined
    if remote is None:
        # then local must be defined
        if local is None:
            raise ValueError("Missing `remote` or `local` arguments")
        mapping = get_local_transcript_gene_mapping(local)
    # remote defined
    else:
        # then local can not be defined
        if local is not None:
            raise ValueError("Only supply `remote` or `local` argument, "
                    "not both.")
        # remote must have 'db'
        if 'db' not in remote:
            raise ValueError("Missing remote database name")
        # and 'annot_src'
        if 'annot' not in remote:
            raise ValueError("Missing remote annotation source name")

        db = remote['db']
        annot = remote['annot']
        if annot not in QUERIES.keys():
            raise ValueError("Invalid annotation source "
                    "name: {0}".format(annot))

        mapping = get_ucsc_transcript_gene_mapping(annot, db, cred=CRED)

    # remove output file if it exists
    if os.path.exists(out_gtf):
        os.remove(out_gtf)

    with track.load(in_gtf, readonly=True) as in_track, \
            track.new(out_gtf, format='gtf') as out_track:
        # since GTF has custom fields, need to set the out_track to use
        # in_track's fields
        out_track.fields = in_track.fields
        for chrom in in_track.chromosomes:
            chrom_rec = in_track.read(chrom)
            out_track.write(chrom, update_gene_id_attr(chrom_rec, mapping))


if __name__ == '__main__':

    usage = __doc__.split('\n\n\n')
    parser = argparse.ArgumentParser(
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=usage[0], epilog=usage[1])

    parser.add_argument('input', type=str, help='Path to input GTF file')
    parser.add_argument('output', type=str, help='Path to output GTF file')
    parser.add_argument('--local', type=str, dest='local', default=None,
            help='Path to transcript ID - gene ID mapping file')
    parser.add_argument('--db', type=str, dest='db', default=None,
            help='UCSC database name to use')
    parser.add_argument('--annot', type=str, dest='annot', default=None,
            choices=QUERIES.keys(), help='UCSC annotation source')
    parser.add_argument('--version', action='version', version='%(prog)s ' +
            __version__)

    args = parser.parse_args()

    remote = None
    if args.db is not None or args.annot is not None:
        remote = {'db': args.db, 'annot': args.annot}

    ucsc_geneid_fix(args.input, args.output, remote=remote, local=args.local)
