#!/bin/bash

## For example: $0 /tmp/inputdownloading /tmp/preprocessing /tmp/processing

## Checking args
if [ $# -ne 3 ]
then
  echo "Usage: $0 inputdownloading_path preprocessing_path processing_path"
  exit 1
fi

## args
INPUTDOWNLOADING_DIR_PATH=$1
PROCESSING_DIR_PATH=$2
PREPROCESSING_DIR_PATH=$3

METADATA_FILE_PATH=$PROCESSING_DIR_PATH/application.metadata
rm -rf $METADATA_FILE_PATH
touch $METADATA_FILE_PATH

MTL_INPUTDOWNLOADING_FILE_PATH=$(find $INPUTDOWNLOADING_DIR_PATH -iname "*_MTL.txt")
STATION_INPUTDOWNLOADING_FILE_PATH=$(find $INPUTDOWNLOADING_DIR_PATH -iname "*_station.csv")

RASTER_ELEVATION_PREPROCESSING_FILE_PATH=$(find $PREPROCESSING_DIR_PATH -iname "elevation.tif")
LAI_PREPROCESSING_FILE_PATH=$(find $PREPROCESSING_DIR_PATH -iname "*_LAI.nc")
SAVI_PREPROCESSING_FILE_PATH=$(find $PREPROCESSING_DIR_PATH -iname "*_SAVI.nc")
NDVI_PREPROCESSING_FILE_PATH=$(find $PREPROCESSING_DIR_PATH -iname "*_NDVI.nc")
LSA_PREPROCESSING_FILE_PATH=$(find $PREPROCESSING_DIR_PATH -iname "*_alb.nc")
LST_PREPROCESSING_FILE_PATH=$(find $PREPROCESSING_DIR_PATH -iname "*_TS.nc")

EF_PROCESSING_FILE_PATH=$(find $PROCESSING_DIR_PATH -iname "*_EF.nc")
ET24H_PROCESSING_FILE_PATH=$(find $PROCESSING_DIR_PATH -iname "*_ET24h.nc")

CURRENT_DATE=$(date)

echo "# Processing (ufcg-sebal) Implementation Metadata" >> $METADATA_FILE_PATH
echo "$CURRENT_DATE # Date" >> $METADATA_FILE_PATH

echo "INPUT" >> $METADATA_FILE_PATH
echo "  INPUTDOWNLOAD" >> $METADATA_FILE_PATH
echo "$MTL_INPUTDOWNLOADING_FILE_PATH # MTL from image" >> $METADATA_FILE_PATH
echo "$STATION_INPUTDOWNLOADING_FILE_PATH # Station data from image" >> $METADATA_FILE_PATH

echo "  PREPROCESSING" >> $METADATA_FILE_PATH
echo "$RASTER_ELEVATION_PREPROCESSING_FILE_PATH # Preprocessed Raster elevation data file path" >> $METADATA_FILE_PATH
echo "$LAI_PREPROCESSING_FILE_PATH # Leaf Area Index data file path" >> $METADATA_FILE_PATH
echo "$SAVI_PREPROCESSING_FILE_PATH # Soil Adjusted Vegetation Index data file path" >> $METADATA_FILE_PATH
echo "$NDVI_PREPROCESSING_FILE_PATH # Normalized Different Vegetation Index data file path" >> $METADATA_FILE_PATH
echo "$LSA_PREPROCESSING_FILE_PATH # Land Surface Albedo data file path" >> $METADATA_FILE_PATH
echo "$LST_PREPROCESSING_FILE_PATH # Land Surface Temperature data file path" >> $METADATA_FILE_PATH

echo "OUTPUT" >> $METADATA_FILE_PATH
echo "$EF_PROCESSING_FILE_PATH # Evapotranspirative fraction data file path" >> $METADATA_FILE_PATH
echo "$ET24H_PROCESSING_FILE_PATH # Evapotranspirative data file path" >> $METADATA_FILE_PATH
