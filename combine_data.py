import sys

import pandas as pd

# IMPORTANT!
# tsv = 0-index
# VCF = 1-indexed

# CONSTANTS: args
# 1st = rna_tsv
# 2nd = sample,
# 3rd = DNA file
# 4th = chromosomal contig

CHUNKSIZE = 10 ** 10
RNA_FILE = sys.argv[0]
SAMPLE = sys.argv[1]
DNA_FILE = sys.argv[2]
CHR = sys.argv[3]
# SAMPLE = 'KEN1092'  # Decide sample name
# RNA_FILE = \
#    '/data/CARDPB/data/NABEC/RNAseq/star_aligned_sorted_bams/' + 
#    'pileups/base_contents/KEN1092fctx/KEN1092fctx_chr21.tsv'
# DNA_FILE = \
#    '/data/CARD_AA/users/johnsonnicl/gene_edits/genotypes/' + \
#    'nabec_chr21_genotypes.tsv'

rna_file = pd.read_csv(RNA_FILE, comment='#', sep='\t', index_col=0)

# First parameter is the replacement, second parameter is your input string
# Memory kinder way to import
# REmoves hyphens in sample name. Be careful about that
def valid(chunks, sample):
    for chunk in chunks:
        mask = chunk['SAMPLE'].str.replace('-', '') == sample
        if mask.all():
            yield chunk
        else:
            yield chunk.loc[mask]
            break


dna_vcf_chunks = pd.read_csv(DNA_FILE,
                             chunksize=CHUNKSIZE,
                             sep='\t',
                             header=0,
                             usecols=['#CHROM', 'POS', 'ID', 'REF', 'ALT',
                                      'QUAL', 'FILTER', 'SAMPLE', 'GT'])
vcf = pd.concat(valid(dna_vcf_chunks, SAMPLE))

# If there are no matches after we already decided on sample overlap then there
# is something wrong
assert len(vcf) >= 1

# testing samples
vcf_t = vcf

vcf_t[['strand1', 'strand2']] = vcf_t['GT'].str.split(pat='/', expand=True)
vcf_t['heterozyg'] = vcf_t['strand1'] != vcf_t['strand2']
vcf_t.drop(columns='GT', inplace=True)

# Worry about indexing
vcf_t.set_index('POS', inplace=True)

# Remove duplicates
drop_index = vcf_t.index.duplicated(keep='first')
vcf_t = vcf_t[~drop_index]
vcf_t.index -= 1

# Join by position (1 has been subtracted from VCF)
# vcf_t = vcf_t.join(rna_file,how='inner',validate='one_to_one')
vcf_t = vcf_t.join(rna_file, how='inner', validate='one_to_one')
vcf_t['COVERAGE'] = vcf_t.loc[:, ['A', 'C', 'T', 'G']].sum(axis=1)
vcf_t = vcf_t[vcf_t['COVERAGE'] >= 30]
vcf_t['alleles'] = vcf_t['strand1'] + vcf_t['strand2']

# Create a pd to output
output_df = pd.DataFrame(columns=['CHR', 'het_or_hom', 'N_in', 'N_out'],
                         index=vcf_t.index)

g_mask = vcf_t['alleles'].str.contains('G')
a_mask = vcf_t['alleles'].str.contains('A')
c_mask = vcf_t['alleles'].str.contains('C')
t_mask = vcf_t['alleles'].str.contains('T')

output_df['N_in'] = 0
output_df.loc[g_mask, 'N_in'] = vcf_t.loc[g_mask, 'G']
output_df.loc[a_mask, 'N_in'] += vcf_t.loc[a_mask, 'A']
output_df.loc[c_mask, 'N_in'] += vcf_t.loc[c_mask, 'C']
output_df.loc[t_mask, 'N_in'] += vcf_t.loc[t_mask, 'T']

g_mask = ~vcf_t['alleles'].str.contains('G')
a_mask = ~vcf_t['alleles'].str.contains('A')
c_mask = ~vcf_t['alleles'].str.contains('C')
t_mask = ~vcf_t['alleles'].str.contains('T')

output_df['N_out'] = 0
output_df.loc[g_mask, 'N_out'] = vcf_t.loc[g_mask, 'G']
output_df.loc[a_mask, 'N_out'] += vcf_t.loc[a_mask, 'A']
output_df.loc[c_mask, 'N_out'] += vcf_t.loc[c_mask, 'C']
output_df.loc[t_mask, 'N_out'] += vcf_t.loc[t_mask, 'T']

output_df['CHR'] = vcf_t['#CHROM']
output_df['het_or_hom'] = vcf_t['heterozyg']

# Change index back to standar vcf 1-index
output_df.index += 1
output_df.reset_index(names='POS', inplace=True)

output_df.to_csv(f"{SAMPLE}_{CHR}_fidelity_counts.tsv",
                 sep='\t', index=False)
