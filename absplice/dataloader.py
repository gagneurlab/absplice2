from collections import defaultdict
import itertools
from typing import List
import pyranges as pr
import numpy as np
import pandas as pd
from tqdm import tqdm
from kipoi.data import SampleIterator
from splicemap.splice_map import SpliceMap
from absplice.utils import read_pangolin, get_pr_coords_pangolin, get_pr_coords_splice_site, get_pr_coords_junctions, pangolin_tissue_specific

try:
    from mmsplice.junction_dataloader import JunctionPSI5VCFDataloader, \
        JunctionPSI3VCFDataloader
    from mmsplice.utils import encodeDNA
except ImportError:
    pass


class SpliceMapMixin:
    def __init__(self, splicemap5=None, splicemap3=None, progress=True):
        self.progress = progress

        if splicemap5 is None and splicemap3 is None:
            raise ValueError(
                '`ref_tables5` and `ref_tables3` cannot be both empty')

        if splicemap5 is not None:
            self.splicemaps5 = self._read_splicemap(splicemap5)
            self.combined_splicemap5 = self._combine_junctions(
                self.splicemaps5)
            self.metadata_splicemap5 = self._splicemap_metadata(self.splicemaps5)
        else:
            self.combined_splicemap5 = None

        if splicemap3 is not None:
            self.splicemaps3 = self._read_splicemap(splicemap3)
            self.combined_splicemap3 = self._combine_junctions(
                self.splicemaps3)
            self.metadata_splicemap3 = self._splicemap_metadata(self.splicemaps3)
        else:
            self.combined_splicemap3 = None

    def _splicemap_metadata(self, splicemaps):
        metadata = defaultdict(list)
        
        cols = ['junction', 'gene_id', 'tissue', 'ref_psi', 'median_n', 'splice_site']
            
        for splicemap in splicemaps:
            df = splicemap.df.copy()
            df = df.rename(columns={'junctions': 'junction'})
            df['tissue'] = splicemap.name 
            itertuples = df[cols].itertuples(index=False)
            
            if self.progress:
                itertuples = tqdm(itertuples)
            
            for row in itertuples:
                metadata[row.junction].append(tuple(row))       
        
        return dict(metadata)
    
    @staticmethod
    def _combine_junctions(splicemaps: List[SpliceMap]):
        columns = ['junctions', 'Chromosome', 'Start', 'End', 'Strand']
        df = pd.concat(
            [s.df[columns] for s in splicemaps]
        ).drop_duplicates(subset='junctions').set_index('junctions')
        return df

    @staticmethod
    def _read_splicemap(path):
        if type(path) is str:
            return [SpliceMap.read_csv(path)]
        elif type(path) is SpliceMap:
            return [path]
        elif type(path) is list:
            return [SpliceMapMixin._read_splicemap(i)[0] for i in path]
        else:
            print(type(path))
            raise ValueError(
                '`splicemap5` or `splicemap3` arguments should'
                ' be list of path to splicemap files'
                ' or `SpliceMap` object')


# class SpliceOutlierDataloader(SpliceMapMixin, SampleIterator):
class SpliceOutlierDataloader(SampleIterator):

    def __init__(self, fasta_file, vcf_file, splicemap_mixin=None):

        import mmsplice
        self.fasta_file = fasta_file
        self.vcf_file = vcf_file
        self.combined_splicemap5 = splicemap_mixin.combined_splicemap5
        self.combined_splicemap3 = splicemap_mixin.combined_splicemap3
        self.splicemaps5 = splicemap_mixin.splicemaps5
        self.splicemaps3 = splicemap_mixin.splicemaps3
        self.metadata_splicemap5 = splicemap_mixin.metadata_splicemap5
        self.metadata_splicemap3 = splicemap_mixin.metadata_splicemap3
        self._generator = iter([])

        if self.combined_splicemap5 is not None:
            self.dl5 = JunctionPSI5VCFDataloader(
                self.combined_splicemap5, fasta_file, vcf_file, encode=False)
            self._generator = itertools.chain(
                self._generator,
                self._iter_dl(self.dl5, self.combined_splicemap5, event_type='psi5'))

        if self.combined_splicemap3 is not None:
            self.dl3 = JunctionPSI3VCFDataloader(
                self.combined_splicemap3, fasta_file, vcf_file, encode=False)
            self._generator = itertools.chain(
                self._generator,
                self._iter_dl(self.dl3, self.combined_splicemap3, event_type='psi3'))

    def _iter_dl(self, dl, intron_annotations, event_type):
        for row in dl:
            junction_id = row['metadata']['exon']['junction']
            ref_row = intron_annotations.loc[junction_id]
            row['metadata']['junction'] = dict()
            row['metadata']['junction']['junction'] = ref_row.name
            row['metadata']['junction']['event_type'] = event_type
            row['metadata']['junction'].update(ref_row.to_dict())
            yield row

    def __next__(self):
        return next(self._generator)

    def __iter__(self):
        return self

    def batch_iter(self, batch_size=32, **kwargs):
        for batch in super().batch_iter(batch_size, **kwargs):
            batch['inputs']['seq'] = self._encode_batch_seq(
                batch['inputs']['seq'])
            batch['inputs']['mut_seq'] = self._encode_batch_seq(
                batch['inputs']['mut_seq'])
            yield batch

    def _encode_batch_seq(self, batch):
        return {k: encodeDNA(v.tolist()) for k, v in batch.items()}
    
class PangolinSpliceMap:
    def __init__(self, df_pangolin, splicemap5=None, splicemap3=None, pyranges_slack=2):
        self.df_pangolin_gain, self.df_pangolin_loss = PangolinSpliceMap._read_pangolin(df_pangolin)
        self.slack = pyranges_slack
        
        if splicemap5 is None and splicemap3 is None:
            raise ValueError("Both PSI5 and PSI3 SpliceMaps cannot be missing")
            
        if splicemap5 is not None:
            df_splicemap5 = PangolinSpliceMap._read_splicemap(splicemap5, 'psi5')
        else:
            df_splicemap5 = None
        if splicemap3 is not None:
            df_splicemap3 = PangolinSpliceMap._read_splicemap(splicemap3, 'psi3')
        else:
            df_splicemap3 = None

        self.df_splicemap = pd.concat([df_splicemap5, df_splicemap3])
        self.join_pangolin_splicemap()
        
    @staticmethod
    def _read_pangolin(path):
        df_pangolin = read_pangolin(path)
        
        df_pangolin_gain = df_pangolin.copy()
        df_pangolin_loss = df_pangolin.copy()
        
        # get coords for the predcited splice sites
        df_pangolin_gain['Chromosome'], df_pangolin_gain['Start'], df_pangolin_gain['End'] = zip(*df_pangolin_gain.apply(
            lambda row: get_pr_coords_pangolin(row, 'gain'), axis=1
        ))
        df_pangolin_loss['Chromosome'], df_pangolin_loss['Start'], df_pangolin_loss['End'] = zip(*df_pangolin_loss.apply(
            lambda row: get_pr_coords_pangolin(row, 'loss'), axis=1
        ))

        return df_pangolin_gain, df_pangolin_loss

    @staticmethod
    def _read_splicemap(path, event_type):
        if type(path) is str:
            return PangolinSpliceMap._read_process_splicemap_file(path, event_type)
        elif type(path) is list:
            return pd.concat([PangolinSpliceMap._read_process_splicemap_file(filename, event_type) for filename in path])
        else:
            raise ValueError(
                f'{event_type} splicemap argument should be list of path to SpliceMap files, not {type(path)}'
            )

    @staticmethod
    def _read_process_splicemap_file(filename, event_type):
        splicemap_cols = [
            'junctions', 
            'gene_id', 
            'splice_site', 
            'ref_psi', 
            'median_n'
        ]

        if type(filename) is SpliceMap:
            sm = filename
        elif type(filename) is str:
            sm = SpliceMap.read_csv(filename)
        else:
            raise ValueError(f"Something went wrong with the SpliceMap type, the passed type is: {type(filename)}")
        df_sm = sm.df[splicemap_cols].copy()
        df_sm['tissue'] = sm.name
        df_sm['event_type'] = event_type
        if df_sm.shape[0] > 0:
            # _df_sm1 = df_sm.copy()
            # _df_sm2 = df_sm.copy()
            
            # _df_sm1['Chromosome'], _df_sm1['Start'], _df_sm1['End'] = zip(*_df_sm1['junctions'].apply(
            #     lambda x: get_pr_coords_junctions(x, 'j1')
            # ))
            # _df_sm2['Chromosome'], _df_sm2['Start'], _df_sm2['End'] = zip(*_df_sm2['junctions'].apply(
            #     lambda x: get_pr_coords_junctions(x, 'j2')
            # ))
            
            # return pd.concat([
            #     _df_sm1,
            #     _df_sm2
            # ])
            df_sm['Chromosome'], df_sm['Start'], df_sm['End'] = zip(*df_sm['splice_site'].apply(
                lambda x: get_pr_coords_splice_site(x)
            ))
        return df_sm

    def join_pangolin_splicemap(self):
        splicemap_cols = [
            'Chromosome', 
            'Start', 
            'End', 
            'junctions', 
            'gene_id', 
            'splice_site', 
            'ref_psi', 
            'median_n', 
            'tissue', 
            'event_type'
        ]

        df_pangolin_loss_joined = pr.PyRanges(self.df_pangolin_loss).join(
            pr.PyRanges(self.df_splicemap[splicemap_cols]), how='left', slack=self.slack
        ).df.drop(columns=['Start_b', 'End_b'])
        # to not lose the pangolin scores when no splice sites overlap
        if df_pangolin_loss_joined['gene_id_b'].unique().tolist() == ['-1']:
            df_pangolin_loss_joined['gene_id_b'] = df_pangolin_loss_joined['gene_id']

        df_pangolin_gain_joined = pr.PyRanges(self.df_pangolin_gain).join(
            pr.PyRanges(self.df_splicemap[splicemap_cols]), how='left', slack=self.slack
        ).df.drop(columns=['Start_b', 'End_b'])
        # to not lose the pangolin scores when no splice sites overlap
        if df_pangolin_gain_joined['gene_id_b'].unique().tolist() == ['-1']:
            df_pangolin_gain_joined['gene_id_b'] = df_pangolin_gain_joined['gene_id']

        # Get the correct gene ids
        df_pangolin_loss_joined = df_pangolin_loss_joined[
            (df_pangolin_loss_joined['gene_id']==df_pangolin_loss_joined['gene_id_b'])
        ].drop(columns='gene_id_b')
        
        df_pangolin_gain_joined = df_pangolin_gain_joined[
            (df_pangolin_gain_joined['gene_id']==df_pangolin_gain_joined['gene_id_b'])
        ].drop(columns='gene_id_b')

        # Replace -1s with NaNs
        str_cols = ['junctions', 'splice_site', 'tissue', 'event_type']
        num_cols = ['ref_psi', 'median_n']

        df_pangolin_loss_joined[str_cols] = df_pangolin_loss_joined[str_cols].replace('-1', None)
        df_pangolin_loss_joined[num_cols] = df_pangolin_loss_joined[num_cols].replace(-1, np.nan)
        
        df_pangolin_gain_joined[str_cols] = df_pangolin_gain_joined[str_cols].replace('-1', None)
        df_pangolin_gain_joined[num_cols] = df_pangolin_gain_joined[num_cols].replace(-1, np.nan)

        # Drop columns
        df_pangolin_gain_joined = df_pangolin_gain_joined.drop(columns=['Chromosome', 'Start', 'End'])
        df_pangolin_loss_joined = df_pangolin_loss_joined.drop(columns=['Chromosome', 'Start', 'End'])

        # Join back together loss and gain dfs
        input_index = [
            'variant',
            'gene_id',
            'tissue',
            'gain_score', 
            'gain_pos', 
            'loss_score', 
            'loss_pos'
        ]

        if 'junctions' in df_pangolin_gain_joined.columns:
            df_pangolin_gain_joined = df_pangolin_gain_joined.rename(columns={'junctions': 'junction'})
        if 'junctions' in df_pangolin_loss_joined.columns:
            df_pangolin_loss_joined = df_pangolin_loss_joined.rename(columns={'junctions': 'junction'})

        df_pangolin_joined = df_pangolin_gain_joined.set_index(input_index).join(
            df_pangolin_loss_joined.set_index(input_index), how='outer', lsuffix='_gain', rsuffix='_loss'
        ).reset_index()

        # Add Pangolin tissue-specific score
        df_pangolin_joined[['pangolin_tissue_score', 'ref_psi_pangolin', 'median_n_pangolin', 'splice_site_pangolin', 'junction_pangolin']] = df_pangolin_joined.apply(
            pangolin_tissue_specific, axis=1, result_type='expand'
        )

        # Generate all unique tissue names
        all_tissues = self.df_splicemap['tissue'].unique()

        # Generate all unique variants - gene IDs
        all_variants = self.df_pangolin_gain.drop(columns=['Chromosome', 'Start', 'End']).drop_duplicates()
        
        # Create all variant-geneID-tissue combinations
        df_variants_tissues = all_variants.merge(pd.DataFrame({'tissue': all_tissues}), how='cross')

        join_index = [
            'variant', 
            'gene_id', 
            'tissue',
            'gain_score', 
            'gain_pos', 
            'loss_score', 
            'loss_pos'
        ]

        self.df_pangolin_splicemap = df_variants_tissues.set_index(join_index).join(
            df_pangolin_joined.set_index(join_index)).reset_index()
