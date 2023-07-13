import pandas as pd

# IMPORTANT!
# tsv = 0-index
# VCF = 1-indexed

file = pd.read_csv('/data/CARDPB/data/NABEC/RNAseq/star_aligned_sorted_bams/pileups/base_contents/UM1792fctx/UM1792fctx_chr21.tsv', comment='#', sep='\t', low_memory=True)

vcf = pd.read_csv('/data/CARD_AA/users/johnsonnicl/gene_edits/genotypes/nabec_chr21_genotypes.tsv', comment='#', sep='\t', low_memory=True)

# testing samples
file_t = file[0:6005000]
vcf_t = vcf[0:100]
