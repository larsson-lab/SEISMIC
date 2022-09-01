#!/bin/bash

# Paths
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
demo_config="$SCRIPT_DIR"'/demo_config.yaml'

# Run SEISMIC on demo data (TCGA UCEC)
echo 'Running SEISMIC on the demo data.'
cd "$SCRIPT_DIR"'/..'
./SEISMIC.R "$demo_config"
