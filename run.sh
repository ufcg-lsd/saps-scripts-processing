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

exit 0

## Exit code
# exit code `0` indicates a successful execution. Any other number indicates failure.
