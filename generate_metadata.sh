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
METADATA_DIR_PATH=$4

METADATA_FILE_PATH=$METADATA_DIR_PATH/metadata.txt
rm -rf $METADATA_FILE_PATH
touch $METADATA_FILE_PATH

MTL_INPUT_FILE_PATH=$(find $INPUT_DIR_PATH -iname "*_MTL.txt")
STATION_INPUT_FILE_PATH=$(find $INPUT_DIR_PATH -iname "*_station.csv")

RASTER_ELEVATION_PREPROCESSING_FILE_PATH=$(find $PREPROCESSING_DIR_PATH -iname "elevation.tif")
LAI_PREPROCESSING_FILE_PATH=$(find $PREPROCESSING_DIR_PATH -iname "*_LAI.nc")
SAVI_PREPROCESSING_FILE_PATH=$(find $PREPROCESSING_DIR_PATH -iname "*_SAVI.nc")
NDVI_PREPROCESSING_FILE_PATH=$(find $PREPROCESSING_DIR_PATH -iname "*_NDVI.nc")
LSA_PREPROCESSING_FILE_PATH=$(find $PREPROCESSING_DIR_PATH -iname "*_alb.nc")
LST_PREPROCESSING_FILE_PATH=$(find $PREPROCESSING_DIR_PATH -iname "*_TS.nc")

EF_OUTPUT_FILE_PATH=$(find $OUTPUT_DIR_PATH -iname "*_EF.nc")
ET24H_OUTPUT_FILE_PATH=$(find $OUTPUT_DIR_PATH -iname "*_ET24h.nc")

CURRENT_DATE=$(date)

echo "# Processing (ufcg-sebal) Implementation Metadata" >> $METADATA_FILE_PATH
echo "$CURRENT_DATE # Date" >> $METADATA_FILE_PATH

echo "INPUT" >> $METADATA_FILE_PATH
echo "  INPUTDOWNLOAD" >> $METADATA_FILE_PATH
echo "$MTL_INPUT_FILE_PATH # MTL from image" >> $METADATA_FILE_PATH
echo "$STATION_INPUT_FILE_PATH # Station data from image" >> $METADATA_FILE_PATH

echo "  PREPROCESSING" >> $METADATA_FILE_PATH
echo "$RASTER_ELEVATION_PREPROCESSING_FILE_PATH # Preprocessed Raster elevation data file path" >> $METADATA_FILE_PATH
echo "$LAI_PREPROCESSING_FILE_PATH # Leaf Area Index data file path" >> $METADATA_FILE_PATH
echo "$SAVI_PREPROCESSING_FILE_PATH # Soil Adjusted Vegetation Index data file path" >> $METADATA_FILE_PATH
echo "$NDVI_PREPROCESSING_FILE_PATH # Normalized Different Vegetation Index data file path" >> $METADATA_FILE_PATH
echo "$LSA_PREPROCESSING_FILE_PATH # Land Surface Albedo data file path" >> $METADATA_FILE_PATH
echo "$LST_PREPROCESSING_FILE_PATH # Land Surface Temperature data file path" >> $METADATA_FILE_PATH

echo "OUTPUT" >> $METADATA_FILE_PATH
echo "$EF_OUTPUT_FILE_PATH # Evapotranspirative fraction data file path" >> $METADATA_FILE_PATH
echo "$ET24H_OUTPUT_FILE_PATH # Evapotranspirative data file path" >> $METADATA_FILE_PATH
