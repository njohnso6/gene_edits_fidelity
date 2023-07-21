import sys
import os
import logging
from pathlib import Path

import numpy as np
import pandas as pd

# IMPORTANT!
# tsv = 0-index
# VCF = 1-indexed

# CONSTANTS: args
# 1st = rna_tsv
# 2nd = sample,
# 3rd = DNA file
# 4th = chromosomal contig
CHUNKSIZE = 10 ** 7
CWD = os.getcwd()

# Test whether using snakemake
isSnakemake = True
try:
    snakemake.input[0]
except NameError:
    isSnakemake = False

if len(sys.argv) > 1:  # If there were command-line arguments
    RNA_FILE = sys.argv[1]
    SAMPLE = sys.argv[2]
    DNA_FILE = sys.argv[3]
    CHR = sys.argv[4]
    LOG = sys.argv[5]
elif isSnakemake:  # if any snakemake arguments
    # Logging

    RNA_FILE = snakemake.input.rna
    SAMPLE = snakemake.wildcards.sample
    DNA_FILE = snakemake.input.dna
    CHR = snakemake.wildcards.chrom
    LOG = snakemake.log[0]
else:
    SAMPLE = 'KEN1092'  # Decide sample name
    RNA_FILE = \
        '/data/CARDPB/data/NABEC/RNAseq/star_aligned_sorted_bams/' + \
        'pileups/base_contents/KEN1092fctx/KEN1092fctx_chr21.tsv'
    DNA_FILE = \
        '/data/CARD_AA/users/johnsonnicl/gene_edits/genotypes/' + \
        'nabec_chr21_genotypes.tsv'
    CHR = 'chr21'
    LOG = "logs/testlog.log"

logging.basicConfig(filename=LOG,
                    format='%(asctime)s %(message)s',
                    filemode='w',
                    level=logging.DEBUG)
logger = logging.getLogger()

logger.debug(msg=CHR)
logger.debug(msg=CWD)
logger.debug(msg=DNA_FILE)
logger.debug(msg=SAMPLE)
logger.debug(msg=RNA_FILE)

# Read in ints as unsigned int16 range 0-65,535
# could probably step down to a single byte if necessary
rna_file = pd.read_csv(RNA_FILE, comment='#', sep='\t',
                       index_col=0, dtype={'A': np.ushort, 'C': np.ushort, 'G' : np.ushort, 'T' : np.ushort})
logger.debug(msg="Read RNA file")

# First parameter is the replacement, second parameter is your input string
# Memory kinder way to import
# REmoves hyphens in sample name. Be careful about that
def valid(chunks, sample):
    for chunk in chunks:
        mask = chunk['SAMPLE'].str.replace('-', '') == sample
        yield chunk.loc[mask]


dna_vcf_chunks = pd.read_csv(DNA_FILE,
                             chunksize=CHUNKSIZE,
                             sep='\t',
                             header=0,
                             usecols=['#CHROM', 'POS', 'ID', 'REF', 'ALT',
                                      'SAMPLE', 'GT'],
                             dtype={'#CHROM': 'category', 'POS': np.uint32,
                                    'REF': 'category', 'ALT': 'category',
                                    'SAMPLE': 'category', 'GT': 'category'})

vcf_t = pd.concat(valid(dna_vcf_chunks, SAMPLE))


# testing samples
# vcf_t = vcf

logger.debug(msg="Imported file")

# If there are no matches after we already decided on sample overlap then there
# is something wrong or the sample and vcf do not overlap
if len(vcf_t) < 1:
    logger.debug(msg="IMPORTANT: No vcf/rna overlap")

    path = f"{CWD}/fidelity_counts/{SAMPLE}"
    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path)

    Path(f"{path}/{SAMPLE}_{CHR}_fidelity_counts.tsv").touch()
    
    with open('no_overlap_files.txt', 'a', encoding='utf8') as fd:
        fd.write(f'\n{SAMPLE}_{CHR}')

    # Append to list of empty files
    with open('empty_files_list.txt', 'a', encoding='utf8') as fd:
        fd.write(f'\n{SAMPLE}_{CHR}')

    # End program
    sys.exit(0)

vcf_t[['strand1', 'strand2']] = vcf_t['GT'].str.split(pat='/', expand=True)
vcf_t['heterozyg'] = vcf_t['strand1'] != vcf_t['strand2']
vcf_t.drop(columns='GT', inplace=True)
logger.debug(msg="Added heterozyg column")

# Worry about indexing
vcf_t.set_index('POS', inplace=True)

# Remove duplicates
drop_index = vcf_t.index.duplicated(keep='first')
vcf_t = vcf_t[~drop_index]
vcf_t.index -= 1
logger.debug(msg="remove duplicates")

# Join by position (1 has been subtracted from VCF)
# vcf_t = vcf_t.join(rna_file,how='inner',validate='one_to_one')
vcf_t = vcf_t.join(rna_file, how='inner', validate='one_to_one')

del rna_file  # Just giving ourselves a bit more memory to work with (possibly)

vcf_t['COVERAGE'] = vcf_t.loc[:, ['A', 'C', 'T', 'G']].sum(axis=1)
vcf_t = vcf_t[vcf_t['COVERAGE'] >= 30]
vcf_t['alleles'] = vcf_t['strand1'] + vcf_t['strand2']
logger.debug(msg="Added coverage and alleles")

# Create a pd to output
output_df = pd.DataFrame(columns=['CHR', 'ID', 'het_or_hom', 'N_in', 'N_out'],
                         index=vcf_t.index)

g_mask = vcf_t['alleles'].str.contains('G')
a_mask = vcf_t['alleles'].str.contains('A')
c_mask = vcf_t['alleles'].str.contains('C')
t_mask = vcf_t['alleles'].str.contains('T')
logger.debug(msg="Masking")

output_df['N_in'] = 0
output_df.loc[g_mask, 'N_in'] = vcf_t.loc[g_mask, 'G']
output_df.loc[a_mask, 'N_in'] += vcf_t.loc[a_mask, 'A']
output_df.loc[c_mask, 'N_in'] += vcf_t.loc[c_mask, 'C']
output_df.loc[t_mask, 'N_in'] += vcf_t.loc[t_mask, 'T']
logger.debug(msg="Added N_in")

g_mask = ~vcf_t['alleles'].str.contains('G')
a_mask = ~vcf_t['alleles'].str.contains('A')
c_mask = ~vcf_t['alleles'].str.contains('C')
t_mask = ~vcf_t['alleles'].str.contains('T')
logger.debug(msg="MAsking")

output_df['N_out'] = 0
output_df.loc[g_mask, 'N_out'] = vcf_t.loc[g_mask, 'G']
output_df.loc[a_mask, 'N_out'] += vcf_t.loc[a_mask, 'A']
output_df.loc[c_mask, 'N_out'] += vcf_t.loc[c_mask, 'C']
output_df.loc[t_mask, 'N_out'] += vcf_t.loc[t_mask, 'T']
logger.debug(msg="Added N_out")

output_df['CHR'] = vcf_t['#CHROM']
output_df['het_or_hom'] = vcf_t['heterozyg']
output_df['ID'] = vcf_t['ID']

# Change index back to standar vcf 1-index
output_df.index += 1
output_df.reset_index(names='POS', inplace=True)


# If df has rows, output df
# else: output a flag file instead indicating empty rows
logger.debug(msg="Writing")
if len(output_df) >= 1:
    # check if the necessary dir exists
    path = f"{CWD}/fidelity_counts/{SAMPLE}"
    isExist = os.path.exists(path)

    if not isExist:
        os.makedirs(path)

    output_df.to_csv(f"{path}/{SAMPLE}_{CHR}_fidelity_counts.tsv",
                     sep='\t', index=False)

    # Append to list of full files
    with open('full_files_list.txt', 'a', encoding='utf8') as fd:
        fd.write(f'\n{SAMPLE}_{CHR}')
    
else:
    path = f"fidelity_counts/{SAMPLE}"
    isExist = os.path.exists(path)
    if not isExist:
        os.makedirs(path)
    logger.debug(msg="No coverage; no output df")
    Path(f"{path}/{SAMPLE}_{CHR}_fidelity_counts.tsv").touch()

    with open('no_coverage_files_list.txt', 'a', encoding='utf8') as fd:
        fd.write(f'\n{SAMPLE}_{CHR}')

    # Append to list of empty files
    with open('empty_files_list.txt', 'a', encoding='utf9') as fd:
        fd.write(f'\n{SAMPLE}_{CHR}')
