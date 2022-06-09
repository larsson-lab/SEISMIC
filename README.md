<img src="https://github.com/larsson-lab/SEISMIC/raw/main/media/SEISMIC_logo.png" width="300" /> <!-- Change this to a raw.github link after making the repo public to avoid broken links if the README is hosted without the logo file in places other than GitHub. -->

## Selection Evidence Inferred from Sample-specific Modelling In Cancer

SEISMIC identifies genomic regions where mutations are positively selected in cancer genomes, by analysing the distribution of mutations across a tumour cohort.

## System Requirements

SEISMIC can be run on a computer or server running Linux, preferably with a decent amount of RAM (at least 32 GB). R is required (version 4.1.0 tested), as well as the following R packages (tested versions in parentheses):
- `tidyverse` (1.3.1)
- `yaml` (2.3.5)
- `foreach` (1.5.2)
- `doSNOW` (1.0.20)
- `stringi` (1.7.6)
- `inline` (0.3.19)
- `Rcpp` (1.0.8.3)
- `fitdistrplus` (1.1.8)
- `data.table` (1.14.2)
- `plyranges` (1.14.0)
- `BSgenome.Hsapiens.UCSC.hg19` (1.4.3) and/or `BSgenome.Hsapiens.UCSC.hg38` (1.4.4)
- `cowplot` (1.1.1)


## Installation

To use the SEISMIC script, simply clone this repository, and run it according to the usage instructions below.
There is an installation script for the required R packages (`demo/install_R_dependencies.R`) that can be sourced from the R console to install the most recent versions. The installation takes approximately 30 minutes if no packages were installed previously.


## Demo
There is a demo script (demo/run_demo.sh) that uses SEISMIC to analyse TCGA UCEC SNVs in gene CDSs. Before running this script, make sure to install the required R packages (see above). SEISMIC will use the settings in `demo/demo_config.yaml`. If the machine running the demo has ample RAM, multi-threading can be enabled by changing the "cores" value in the config file to speed up analysis. For reference, 5-10 GB of RAM may be used per thread.

The demo script runs SEISMIC, which first annotates mutation effects in gene CDSs (10 min using 16 GB RAM on our system with 16 threads available), and then runs the test (1 h 40 min using \~75 GB RAM running on 10 threads on our system). The expected output is the file `demo/output/SEISMIC_demo_UCEC_WXS_result__unadjusted_bkg_mutrate_cohort_sig_UCEC_min_3_muts_10000_sims.tsv`, with genes listed by P value and the same results as in the UCEC analysis in the paper.


## Usage

SEISMIC is run from the command line, with a yaml config file as input argument.

```bash
Rscript SEISMIC.R config.yaml
```

### Config file

The SEISMIC script requires a yaml config file to run. A template config file with default values is supplied in `config/config_template.yaml`, annotated with explanations of the variables.


### Preparation

Mutations (hg19 or hg38 coordinates) can either be in MAF format, or tab-separated with column names at the top of the file, with the following columns:
- `seqnames`: chromsome name (e.g., chr1)
- `start`: 1-based position of mutation
- `end`: 1-based position of mutation
- `strand`: strand of the refnuc/varnuc 
- `cancer`: cancer type. Multiple cancer types in one file will filter to the specified cancer type.
- `refnuc`: reference nucleotide
- `varnuc`: variant nucleotide
- `sampleID`: patient ID

A file with region definitions for the regions to be analysed is required. The files used for gene CDSs, promoters, and introns in the paper are available in the data directory.

When running SEISMIC with mutation expectation scaling to replication timing or expression, the relevant files are available in the data directory. Please note that this is only relevant when analysing recurrence as well as cohort distribution, and not required for normal use as in the paper.

