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

# fake files
touch $PROCESSING_DIR_PATH/dados.csv
touch $PROCESSING_DIR_PATH/error.log
touch $PROCESSING_DIR_PATH/LC82150652015174LGN00_cpu_usage.txt
touch $PROCESSING_DIR_PATH/LC82150652015174LGN00_disk_usage.txt
touch $PROCESSING_DIR_PATH/LC82150652015174LGN00_EF.nc
touch $PROCESSING_DIR_PATH/LC82150652015174LGN00_ET24h.nc
touch $PROCESSING_DIR_PATH/LC82150652015174LGN00_mem_usage.txt
touch $PROCESSING_DIR_PATH/out.log
touch $PROCESSING_DIR_PATH/stage.metadata
touch $PROCESSING_DIR_PATH/timestamp

exit 0

## Exit code
# exit code `0` indicates a successful execution. Any other number indicates failure.
