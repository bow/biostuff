#!/usr/bin/env sh

# get_ucsc_ensembl_enst_ensg_map.sh
#
# Simple script for downloading Ensembl IDs of canonical UCSC transcripts.
#
# Usage: ./get_ucsc_ensembl_enst_ensg_map.sh > {your_output_file}

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A --column-names=FALSE << QUERY
use hg19;
SELECT ensGene.name, ensGene.name2
    FROM ensGene
    ORDER BY LENGTH(ensGene.chrom), ensGene.chrom, ensGene.txStart;
QUERY
