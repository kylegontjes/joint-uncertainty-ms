#!/usr/bin/env python3

from snakemake.shell import shell
import os 
import sys
import pandas as pd 
from collections import defaultdict
import Bio
from Bio import AlignIO
from Bio import SeqIO  
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

outdir = snakemake.params[0]
prefix = snakemake.params[1]
alignment = snakemake.params[2] 
gff_file = snakemake.input[0]
 
# read in alignment and gubbins gff file
gff = pd.read_csv(gff_file, sep='\t', skiprows=2, header=None)
aln = AlignIO.read(alignment, 'fasta')

# get indices/positions of recombinant regions identified by gubbins
#print('Getting recombinant positions.')
recomb_regions = defaultdict(list)

for row in gff.iterrows():
    start = row[1][3]
    end = row[1][4]
    region = list(range(start, end+1))
    taxa = row[1][8].split(';')[2]
    taxa = taxa.replace('taxa=\"', '')
    taxa = taxa.replace('\"', '')
    taxa = list(taxa.split())
    for isolate in taxa:
        for position in region:
            recomb_regions[isolate].append(position)


# mask indices/positions of recombinant regions identified by gubbins
#print('Masking recombinant positions in whole genome alignment.')
sample_masked_indices = defaultdict(list)
new_aln = list()

for record in aln:
    seq_str = list(str(record.seq))
    masked_indices = recomb_regions.get(record.id, [])
    for index in masked_indices:
        seq_str[index-1] = 'N'
    seq_str = ''.join(seq_str)
    new_record = SeqRecord(Seq(seq_str), id=record.id, description="")
    sample_masked_indices[record.id] = masked_indices
    new_aln.append(new_record)

# write new FASTA file with recombinant regions masked
fasta_outfile = os.path.join(outdir, prefix + "_gubbins_masked.fa")
text_outfile = os.path.join(outdir, prefix + "_masked_recomb_positions.txt")
var_site_outfile = os.path.join(outdir,prefix + "_gubbins_masked_var_sites.fa")

#print('Writing', fasta_outfile)
with open(fasta_outfile, 'w') as handle:
    SeqIO.write(new_aln, handle, 'fasta')

# Write text file with list of recombinant sites for each genome
#print('Writing', text_outfile)
with open(text_outfile, 'w') as handle:
    for sample, positions in sample_masked_indices.items():
        line = str(sample) + '\t' + ','.join(map(str, positions)) + '\n'
        handle.write(line)

# Get variant sites and write to fasta file using snp-sites
#print('Getting variant sites using snp-sites.')
cmd = 'snp-sites ' + fasta_outfile + \
      ' -m ' + ' -o ' + var_site_outfile 

os.system(cmd)
