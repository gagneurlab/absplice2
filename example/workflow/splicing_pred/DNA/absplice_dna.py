import pandas as pd
import numpy as np
import pickle
from tqdm import tqdm
from splicemap.splice_map import SpliceMap
tqdm.pandas()

def get_abs_max_rows(df, groupby, max_col, dropna=True):
    return df.reset_index() \
        .sort_values(by=max_col, key=abs, ascending=False) \
        .drop_duplicates(subset=groupby) \
        .set_index(groupby)

df_pangolin = pd.read_csv(snakemake.input['pangolin_splicemap'])
df_mmsplice_splicemap = pd.read_csv(snakemake.input['mmsplice_splicemap'])

df_mmsplice_splicemap = df_mmsplice_splicemap[[
    'variant', 'gene_id', 'tissue',
    'ref_psi', 'median_n', 'delta_logit_psi', 'delta_psi',
    'junction', 'event_type', 'splice_site'
]]
df_pangolin = df_pangolin[[
    'variant', 'gene_id', 'tissue', 
    'gain_score', 'gain_pos', 'loss_score', 'loss_pos', 
    'pangolin_tissue_score', 'ref_psi_pangolin', 'median_n_pangolin',
]]

join_index = ['variant', 'gene_id', 'tissue']
df_pangolin = df_pangolin.drop_duplicates()
df_mmsplice_splicemap = df_mmsplice_splicemap.drop_duplicates()
df_joined = df_pangolin.set_index(join_index).join(
        df_mmsplice_splicemap.set_index(join_index), how='outer'
    ).reset_index()
    
# clip pangolin gain score, due to feature contribution function (not enough gain score outliers in training data)    
df_joined['gain_score_original'] = df_joined['gain_score'].copy()
df_joined['gain_score'] = df_joined['gain_score'].clip(upper=0.7)

del df_mmsplice_splicemap
del df_pangolin

features = [
    'delta_logit_psi',
    'delta_psi',
    'gain_score',
    'loss_score',
    'median_n',
    'median_n_pangolin',
]

absplice_model = '../../absplice/precomputed/AbSplice_2_DNA.pkl'

model = pickle.load(open(absplice_model, 'rb'))
df_features = df_joined[features].fillna(0).copy()
df_joined['AbSplice_DNA'] = model.predict_proba(df_features)[:, 1]

df_joined = get_abs_max_rows(df_joined, ['variant', 'gene_id', 'tissue'], 'AbSplice_DNA')

if 'index' in df_joined.columns:
    df_joined = df_joined.drop(columns=['index'])

df_joined['gain_score'] = df_joined['gain_score_original'].copy()
df_joined = df_joined.drop(columns='gain_score_original')

df_joined.reset_index().to_csv(snakemake.output['absplice_dna'], index=False)