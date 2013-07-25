#!/usr/bin/env sh

# get_ucsc_ensembl_bed.sh
#
# Simple script for downloading Ensembl transcripts from UCSC as a BED file.
#
# Usage: ./get_ucsc_ensembl_gene_bed.sh > {your_output_file}

echo track name=\"UCSC Ensembl transcripts\" type=bed description=\"downloaded on `date +%F`\"
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A --column-names=FALSE << QUERY
use hg19;
SELECT chrom, txStart, txEnd, name, score, strand
    FROM ensGene
    ORDER BY LENGTH(chrom), chrom, txStart, name;
QUERY
