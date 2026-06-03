# scan2 dog/canFam3 support patch

This repository contains scripts and notes for running `scan2` on the dog genome assembly `canFam3`. The original `scan2` workflow does not fully support dog/canFam3, so several functions and main scripts were modified to enable somatic mutation signature analysis in dog single-cell whole-genome sequencing data.

## Installation

Create a conda environment and install `scan2`:

```bash
mamba create -n scan2 -c conda-forge -c bioconda -c jluquette -c dranew -c soil scan2
conda activate scan2
```

## Install SigProfilerMatrixGenerator dog reference

Install the dog reference annotation required by `SigProfilerMatrixGenerator`:

```python
from SigProfilerMatrixGenerator import install as genInstall
genInstall.install("dog", rsync=False, bash=True)
```

If installation fails because of firewall or network restrictions, the reference files can be downloaded manually:

```bash
cd [scan2_env_path]/lib/python3.10/site-packages/SigProfilerMatrixGenerator/references/chromosomes/tsb

wget ftp://alexandrovlab-ftp.ucsd.edu/pub/tools/SigProfilerMatrixGenerator/dog.tar.gz
tar -xzvf dog.tar.gz
```

Replace `[scan2_env_path]` with the path to the conda environment containing `scan2`.

## Install dog BSgenome package in R

The patched `scan2` scripts require the dog genome package `BSgenome.Cfamiliaris.UCSC.canFam3`.

In R:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BSgenome.Cfamiliaris.UCSC.canFam3")
```

## Patch scan2 to support canFam3

The modified `scan2` R functions for dog/canFam3 support are:

```bash
scan2_patch_rScripts.R
```

The modified `scan2` main script, which avoids errors caused by unsupported `canFam3` genome labels, is currently located at:

```bash
scan2_patch_mainScripts.py
```

These patched scripts should be used to replace or override the corresponding original `scan2` scripts in the installed conda environment.  
R script replaced to [scan2_env_path]/lib/R/library/scan2/R/scan2
Main script replaced to [scan2_env_path]/bin/scan2

## Analysis and plotting

The R script used for downstream analysis and plotting is provided as:

```bash
analysis_and_plot.R
```

This script includes the parameters used to generate the final figures and summary plots.

## Notes

* This patch was designed specifically for dog/canFam3 analyses.
* The SigProfilerMatrixGenerator reference is installed under the `"dog"` genome label.
* The R-based genome sequence support is provided by `BSgenome.Cfamiliaris.UCSC.canFam3`.
* Users should verify that chromosome naming conventions are consistent across mutation files, reference genome files, and annotation files before running the workflow.

## Citation

If using this patched workflow, please cite the original `scan2` software and the following paper:  
Class & Yu et al., Dog and human neurons share somatic mutation rates and patterns. Under Revision.

