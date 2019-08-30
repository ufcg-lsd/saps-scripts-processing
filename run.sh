#!/bin/bash

## Objetive: Run the algorithm for the input data downloaded and/or preprocessed
## For example: $0 /tmp/input /tmp/output /tmp/preprocessing /tmp/metadata

## Checking args
if [ $# -ne 4 ]
then
  echo "Usage: $0 input_path output_path preprocessing_path metadata_path"
  exit 1
fi

## Capture args
INPUT_DIR_PATH=$1
OUTPUT_DIR_PATH=$2
PREPROCESSING_DIR_PATH=$3
METADATAPATH_DIR_PATH=$4

# Global variables
SANDBOX=$(pwd)
R_EXEC_DIR=$SANDBOX/src
SCRIPTS_DIR=$SANDBOX/scripts
TMP_DIR_PATH=/tmp
MAX_TRIES=2

R_ALGORITHM_PROCESSING_VERSION=sebal.R

echo "Step 1. Capture MTL and station path"

files=($INPUT_DIR_PATH/*_MTL.txt)
for file in "${files[@]}"
do
  filename="${file##*/}"
  filenameWithoutExtension="${filename%_MTL*}"
  IMAGE_NAME="$filenameWithoutExtension"
done

IMAGE_MTL_PATH=$INPUT_DIR_PATH/$IMAGE_NAME"_MTL.txt"
IMAGE_STATION_FILE_PATH=$INPUT_DIR_PATH/$IMAGE_NAME"_station.csv"

echo "MTL path: $IMAGE_MTL_PATH"
echo "Station path: $IMAGE_STATION_FILE_PATH"

echo "Step 2. Creating dados.csv for image $IMAGE_NAME"

echo "File images;MTL;Path Prepoc;File Station Weather;Path Output" > $R_EXEC_DIR/dados.csv
echo "$INPUT_DIR_PATH;$PREPROCESSING_DIR_PATH;$IMAGE_MTL_PATH;$IMAGE_STATION_FILE_PATH;$OUTPUT_DIR_PATH" >> $R_EXEC_DIR/dados.csv

echo "Step 3. Starting CPU, disk and Memory collect..."

bash $SCRIPTS_DIR/collect-cpu-usage.sh $(pidof R) | tee $OUTPUT_DIR_PATH/$IMAGE_NAME"_cpu_usage.txt" > /dev/null &
bash $SCRIPTS_DIR/collect-memory-usage.sh $(pidof R) | tee $OUTPUT_DIR_PATH/$IMAGE_NAME"_mem_usage.txt" > /dev/null &
bash $SCRIPTS_DIR/collect-disk-usage.sh $(pidof R) | tee $OUTPUT_DIR_PATH/$IMAGE_NAME"_disk_usage.txt" > /dev/null &

echo "Step 4. Executing R script"

for i in `seq $MAX_TRIES`
  do
  
  rm -r $TMP_DIR_PATH/* $METADATA_DIR_PATH/*
  bash $SCRIPTS_DIR/executeRScript.sh $R_EXEC_DIR/$R_ALGORITHM_PREPROCESSING_VERSION $R_EXEC_DIR $TMP_DIR_PATH $METADATA_DIR_PATH
  PROCESS_OUTPUT=$?
  
  echo "executeRScript_process_output=$PROCESS_OUTPUT"
  if [ $PROCESS_OUTPUT -eq 0 ]
  then
    echo "Number of tries: $i"
    break
  elif [ $PROCESS_OUTPUT -eq 124 ] && [ $i -ge $MAX_TRIES ]
  then
    exit 124
  else
    if [ $i -ge $MAX_TRIES ]
    then
      echo "Number of tries: $i"
      exit 1
    fi
  fi
done

echo "Step 5. Killing collect CPU, disk and Memory scripts"

ps -ef | grep collect-cpu-usage.sh | grep -v grep | awk '{print $2}' | xargs kill
ps -ef | grep collect-memory-usage.sh | grep -v grep | awk '{print $2}' | xargs kill
ps -ef | grep collect-disk-usage.sh | grep -v grep | awk '{print $2}' | xargs kill

## Exit code
# The `run.sh` script should have the following return pattern:
# - `0` represents a successful execution.
# And any other exit code will be considered failed.
