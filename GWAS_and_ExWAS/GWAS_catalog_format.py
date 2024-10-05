# Follow instructions from https://www.ebi.ac.uk/gwas/docs/submission
# docker run 0b2f2cc1442f --help

import pandas as pd
import numpy as np
import os

inputdir='/Volumes/acasis/cardiac/48K/gwas_aha_t64/gwas_raw_p2/'
outputdir='/Volumes/acasis/cardiac/48K/gwas_aha_t64/gwas_catalog/'

for fpath in [x for x in os.listdir(inputdir) if x.endswith('.regenie.gz')]:
    path=inputdir + fpath
    df=pd.read_csv(path,sep='\t')
    # select LOG10P > -log10(5e-8)
    df = df[df['LOG10P'] > -np.log10(1e-2)]
    print(fpath, df.shape[0])
    # rename columns
    df.rename(columns={
        'CHROM':'chromosome',
        'GENPOS':'base_pair_location',
        'ALLELE0':'other_allele',
        'ALLELE1':'effect_allele',
        'BETA':'beta',
        'SE':'standard_error',
        'A1FREQ':'effect_allele_frequency',
        'LOG10P':'neg_log_10_p_value',
        'ID':'rs_id',
        'INFO':'info',
        'N':'n',
        'TEST':'test',
    },inplace=True)

    # order columns as above
    df=df[['chromosome',
        'base_pair_location',
        'effect_allele',
        'other_allele',
        'beta',
        'standard_error',
        'effect_allele_frequency',
        'neg_log_10_p_value',
        'rs_id',
        'info',
        'n','test']]

    # save  to tsv
    path= outputdir + fpath.replace('.regenie.gz','.gwas_catalog.tsv.gz')
    df.to_csv(path,sep='\t',index=False)