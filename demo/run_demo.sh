#!/bin/bash

# Paths
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
annotate_mut_effect_script="$SCRIPT_DIR"'/../data/annotate_mutation_effects_in_gene_cds.R'
demo_config="$SCRIPT_DIR"'/demo_config.yaml'

# Run SEISMIC on demo data (TCGA UCEC)
echo 'Running SEISMIC on the demo data.'
cd "$SCRIPT_DIR"'/..'
Rscript SEISMIC.R "$demo_config"
