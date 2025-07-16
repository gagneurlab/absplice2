import pandas as pd
import gzip
import pysam

vcf_path = snakemake.input['pangolin_raw']

df = {
    'variant': [],
    'gene_id': [],
    'gain_score': [],
    'gain_pos': [],
    'loss_score': [],
    'loss_pos': [],
}

vcf_in = pysam.VariantFile(vcf_path, "r")

for variant in vcf_in:
    if 'Pangolin' in variant.info:
        chrom = variant.chrom
        pos = int(variant.pos)
        ref = variant.ref
        alt = variant.alts[0]
        assert variant.info['Pangolin'][0].count('ENS') == 1
        per_gene_string_scores = variant.info['Pangolin']
        for gene_scores in per_gene_string_scores:
            df['variant'].append(f'{chrom}:{pos}:{ref}>{alt}')
            gene_id = gene_scores.split('|')[0].split('.')[0]
            df['gene_id'].append(gene_id)
            gain_score = float(gene_scores.split('|')[1].split(':')[1])
            gain_pos = int(gene_scores.split('|')[1].split(':')[0])
            loss_score = float(gene_scores.split('|')[2].split(':')[1])
            loss_pos = int(gene_scores.split('|')[2].split(':')[0])
            df['gain_score'].append(gain_score)
            df['gain_pos'].append(gain_pos)
            df['loss_score'].append(loss_score)
            df['loss_pos'].append(loss_pos)

vcf_in.close()
df = pd.DataFrame(df)

# # make a new column with the abs max of gain and loss scores
# df['Pangolin_max_score'] = df[['gain_score', 'loss_score']].abs().max(axis=1)

df.to_csv(snakemake.output['pangolin_csv'], index=False)