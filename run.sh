#!/bin/bash

## This script calculates the evapotranspiration data based on the previous phases of the SAPS pipeline

## Checking args
if [ $# -ne 4 ]
then
  echo "Usage: $0 /tmp/teste landsat_X PPPRRR YYYY-MM-DD"
  exit 1
fi

## args
ROOT_DIR=$1
IMAGE_DATASET=$2
IMAGE_PATHROW=$3
IMAGE_DATE=$4

# folders
PROCESSING_DIR_PATH=$ROOT_DIR/processing

# repo with results
REPO=http://www2.lsd.ufcg.edu.br/~thiagoyes/saps/nop-download-files/processing

cd $PROCESSING_DIR_PATH

# download files
wget $REPO/dados.csv
wget $REPO/error.log
wget $REPO/out.log
wget $REPO/stage.metadata
wget $REPO/timestamp
wget $REPO/LC82150652015174LGN00_EF.nc
wget $REPO/LC82150652015174LGN00_ET24h.nc
wget $REPO/LC82150652015174LGN00_cpu_usage.txt
wget $REPO/LC82150652015174LGN00_disk_usage.txt
wget $REPO/LC82150652015174LGN00_mem_usage.txt

exit 0

## Exit code
# exit code `0` indicates a successful execution. Any other number indicates failure.
# In particular, `3` indicates that a Landsat image was not found for the given paramenters.
