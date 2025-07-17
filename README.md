# AbSplice 2

## Installation with creating a conda environment

Clone git repo:
```
git clone https://github.com/gagneurlab/absplice2.git
```

cd into repo directory:
```
cd absplice2
```

Install conda environment (**Note**: users who already have the AbSplice 1 environment installed can skip the steps below and use it):
```
# Recommended if you have mamba installed
mamba env create -f environment.yaml
# otherwise
conda env create -f environment.yaml
```
Activate conda environment:
```
conda activate absplice
```
Install modules from absplice:
```
pip install -e .
```

## Example use case
The [example](https://github.com/gagneurlab/absplice2/tree/master/example) folder contains a snakemake workflow to generate AbSplice predictions, given a vcf file and a fasta file (either for hg19 or hg38, will be downloaded automatically). \
The snakemake workflow will download precomputed SpliceMaps from Zenodo and run AbSplice based on these annotations.
To generate predictions run:
```
cd example/workflow
python -m snakemake -j 1 --use-conda
```
### AbSplice-DNA:
To run the workflow on your own data do the following:

- Store all (or provide a symlink to) vcf files for analysis in [`data/resources/analysis_files/vcf_files/`](https://github.com/gagneurlab/absplice2/tree/master/example/data/resources/analysis_files/vcf_files).

- Specify the genome version that you are going to use (hg19 or hg38) in the field `genome` of the [config](https://github.com/gagneurlab/absplice2/blob/master/example/workflow/config.yaml#L4) file.

- In the field `splicemap_tissues` of the [config](https://github.com/gagneurlab/absplice2/blob/master/example/workflow/config.yaml#L21) file you can uncomment the tissues that AbSplice will use to generate predictions.