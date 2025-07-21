# AbSplice 2
AbSplice 2 is a method that predicts aberrant splicing across human tissues and developmental stages, as described in [Wagner, Çelik et al., Nature Genetics 2023](https://www.nature.com/articles/s41588-023-01373-3), [Wagner et al., biorxiv 2025](https://www.biorxiv.org/content/10.1101/2025.07.16.665183v1).
## Installation with creating a conda environment

Clone git repo:
```
git clone https://github.com/gagneurlab/absplice2.git
```

cd into repo directory:
```
cd absplice2
```

> _**Note**: users who already have the AbSplice 1 environment installed can skip the steps below and continue from **Install modules from absplice**_

Install conda environment:
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
If you want to reuse your AbSplice 1 environment, you also need to install the updated [splicemap](https://github.com/gagneurlab/splicemap/tree/dev) package:
```
pip install git+https://github.com/gagneurlab/splicemap.git@dev
```

## Example use case
The [example](https://github.com/gagneurlab/absplice2/tree/main/example) folder contains a snakemake workflow to generate AbSplice predictions, given a vcf file and a fasta file (either for hg19 or hg38, will be downloaded automatically). \
The snakemake workflow will download precomputed SpliceMaps from Zenodo and run AbSplice based on these annotations.
To generate predictions run:
```
cd example/workflow
python -m snakemake -j 1 --use-conda
```
### AbSplice-DNA:
To run the workflow on your own data do the following:

- Store all (or provide a symlink to) vcf files for analysis in [`data/resources/analysis_files/vcf_files/`](https://github.com/gagneurlab/absplice2/tree/main/example/data/resources/analysis_files/vcf_files).

- Specify the genome version that you are going to use (hg19 or hg38) in the field `genome` of the [config](https://github.com/gagneurlab/absplice2/blob/main/example/workflow/config.yaml#L4) file.

- In the field `splicemap_tissues` of the [config](https://github.com/gagneurlab/absplice2/blob/dev/example/workflow/config.yaml#L24) file you can uncomment the tissues that AbSplice will use to generate predictions.

- **NEW!** If you want to generate developmental predictions, set the field [`devAbSplice`](https://github.com/gagneurlab/absplice2/blob/dev/example/workflow/config.yaml#L10) to `True` and uncomment the desired tissues and timepoints in the fields [`dev_splicemap_tissues`](https://github.com/gagneurlab/absplice2/blob/dev/example/workflow/config.yaml#L75) and [`dev_splicemap_timepoints`](https://github.com/gagneurlab/absplice2/blob/dev/example/workflow/config.yaml#L84) of the config file.
