#!/bin/bash

R_ALGORITHM_PATH=$1
R_EXEC_DIR=$2
TMP_DIR_PATH=$3
METADATA_DIR_PATH=$4

TIMEOUT=1200

timeout $TIMEOUT Rscript $R_ALGORITHM_PATH $R_EXEC_DIR $TMP_DIR_PATH > $METADATA_DIR_PATH/out.log 2> $METADATA_DIR_PATH/error.log
PROCESS_OUTPUT=$?

echo "RScript_process_output=$PROCESS_OUTPUT"
if [ $PROCESS_OUTPUT -eq 124 ]
then
  exit 124
elif [ $PROCESS_OUTPUT -ne 0 ]
then
  exit 1
else
  exit 0
fi
