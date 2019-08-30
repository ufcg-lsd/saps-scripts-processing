########################################################################################
#                                                                                      #
#                         EU BRAZIL Cloud Connect                                      #
#                                                                                      #
#                                                                                      #
########################################################################################

options(echo=TRUE)
rm(list=ls())

library(R.utils)
library(raster)
library(rgdal)
library(maptools)
library(ncdf4)
library(sp)
library(snow)
#library(remotes)
library(sebkc)

##### Opening files ######

args <- commandArgs(trailingOnly=TRUE)
WD <- args[1]
setwd(WD) # Working Directory

# Changing raster tmpdir (Precisa?)
rasterOptions(tmpdir=args[2])

dados <- read.csv("dados.csv", sep=";", stringsAsFactors=FALSE)
MTL <- read.table(dados$MTL[1], skip=0, nrows=140, sep="=", quote = "''", as.is=TRUE) # Reading MTL File

fic <- substr(MTL$V2[MTL$V1 == grep(pattern="LANDSAT_SCENE_ID", MTL$V1, value=T)], 3, 23)
Dia.juliano <- as.numeric(substr(fic, 14, 16))	#Julian Day

sun_elevation <- as.numeric(MTL$V2[MTL$V1 == grep(pattern="SUN_ELEVATION", MTL$V1, value=TRUE)])

fic.sw <- dados$File.Station.Weather[1]
table.sw <- (read.csv(fic.sw, sep=";", header=FALSE, stringsAsFactors=FALSE))
hour.image <- (as.numeric(substr(MTL$V2[MTL$V1 == grep(pattern="SCENE_CENTER_TIME", MTL$V1, value=T)], 3, 4))+
                 as.numeric(substr(MTL$V2[MTL$V1 == grep(pattern="SCENE_CENTER_TIME", MTL$V1, value=T)], 6, 7))/60)*100
hour.image.station<-which.min(abs(table.sw$V3[]-hour.image))

fic.preproc <- dados$Path.Prepoc[1]

raster.elevation <- raster(paste(fic.preproc, "/elevation.tif", sep=""))
albedo <- raster(paste(fic.preproc, "/", fic, "_alb.nc", sep=""))
Ts <- raster(paste(fic.preproc, "/", fic, "_TS.nc", sep=""))
NDVI <- raster(paste(fic.preproc, "/", fic, "_NDVI.nc", sep=""))
LAI <- raster(paste(fic.preproc, "/", fic, "_LAI.nc", sep=""))
SAVI <- raster(paste(fic.preproc, "/", fic, "_SAVI.nc", sep=""))

output.path<-paste(dados$Path.Output[1], "/", fic, ".nc", sep="")

# Time image
acquired_date <- as.Date(MTL$V2[MTL$V1==grep(pattern="DATE_ACQUIRED", MTL$V1, value=TRUE)])
daysSince1970 <- as.numeric(acquired_date)
tdim <- ncdim_def("time", "days since 1970-1-1", daysSince1970, unlim=TRUE, create_dimvar=TRUE, "standard", "time")

##### Opening files ######

azom<- -3    #Par?metro para o Zom imagem
bzom<- 6.47  #Par?metro para o Zom imagem
zom<- exp(azom+bzom*NDVI)

series=sebal(albedo, Ts, NDVI, SAVI, welev=385, xyhot = "auto", xycold = "auto",
    DOY = Dia.juliano, sunelev=sun_elevation, zx = 10,
    u = table.sw$V6[hour.image.station], zomw = 0.2, zom = zom,
    LAI = LAI, DEM=raster.elevation, lapse = 0.0065, Rn24 = NULL, ETr = NULL,
    ETr24 = NULL, wmo = NULL, airport = NULL, Krs = 0.19,
    surface = "grass", latitude = table.sw$V4[1], t1 = 1, 
    time = 12.5, Lz = NULL, Lm = NULL, model = "SEBAL", iter.max = 7, clip = NULL,
    folder = NULL)

output.evapo<-stack(series$EF, series$ETins,series$Rn,series$G)
names(output.evapo)<-c('EF','ET24h',"Rn","G")
writeRaster(output.evapo, output.path, overwrite=TRUE, format="CDF", varname=fic, varunit="daily", longname=fic, xname="lon", yname="lat", bylayer=TRUE, suffix="names")


# Opening old G NetCDF
var_output <- paste(dados$Path.Output[1], "/", fic, "_G.nc", sep="")
nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)

# Getting lat and lon values from old NetCDF
oldLat <- ncvar_get(nc, "lat", start=1, count=raster.elevation@nrows)
oldLon <- ncvar_get(nc, "lon", start=1, count=raster.elevation@ncols)

# Defining latitude and longitude dimensions
dimLatDef <- ncdim_def("lat", "degrees", oldLat, unlim=FALSE, longname="latitude")
dimLonDef <- ncdim_def("lon", "degrees", oldLon, unlim=FALSE, longname="longitude")

# New G file name
file_output <- paste(dados$Path.Output[1], "/", fic, "_G.nc", sep="")
oldGValues <- ncvar_get(nc, fic)
newGValues <- ncvar_def("G", "daily", list(dimLonDef, dimLatDef, tdim), longname="G", missval=NaN, prec="double")
nc_close(nc)
newGNCDF4 <- nc_create(file_output, newGValues)
ncvar_put(newGNCDF4, "G", oldGValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
nc_close(newGNCDF4)

proc.time()

# Opening old Rn NetCDF
var_output <- paste(dados$Path.Output[1], "/", fic, "_Rn.nc", sep="")
nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)

# New Rn file name
file_output <- paste(dados$Path.Output[1], "/", fic, "_Rn.nc", sep="")
oldRnValues <- ncvar_get(nc, fic)
newRnValues <- ncvar_def("Rn", "daily",list(dimLonDef, dimLatDef, tdim), longname="Rn", missval=NaN, prec="double")
nc_close(nc)
newRnNCDF4 <- nc_create(file_output, newRnValues)
ncvar_put(newRnNCDF4, "Rn", oldRnValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
nc_close(newRnNCDF4)

proc.time()


# Opening old EF NetCDF
var_output<-paste(dados$Path.Output[1],"/",fic,"_EF.nc",sep="")
nc<-nc_open(var_output, write=TRUE,readunlim=FALSE,verbose=TRUE,auto_GMT=FALSE,suppress_dimvals=FALSE)

# New EF file name
file_output<-paste(dados$Path.Output[1],"/",fic,"_EF.nc",sep="")
oldEFValues<-ncvar_get(nc,fic)
newEFValues<-ncvar_def("EF","daily",list(dimLonDef,dimLatDef,tdim),longname="EF",missval=NaN,prec="double")
nc_close(nc)
newEFNCDF4<-nc_create(file_output,newEFValues)
ncvar_put(newEFNCDF4,"EF",oldEFValues,start=c(1,1,1),count=c(raster.elevation@ncols,raster.elevation@nrows,1))
nc_close(newEFNCDF4)

print(proc.time())

# Opening old ET24h NetCDF
var_output<-paste(dados$Path.Output[1],"/",fic,"_ET24h.nc",sep="")
nc<-nc_open(var_output, write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)

# New ET24h file name
file_output<-paste(dados$Path.Output[1],"/",fic,"_ET24h.nc",sep="")
oldET24hValues<-ncvar_get(nc,fic)
newET24hValues<-ncvar_def("ET24h","daily", list(dimLonDef, dimLatDef, tdim), longname="ET24h", missval=NaN, prec="double")
nc_close(nc)
newET24hNCDF4<-nc_create(file_output,newET24hValues)
ncvar_put(newET24hNCDF4, "ET24h", oldET24hValues, start=c(1, 1, 1), count=c(raster.elevation@ncols, raster.elevation@nrows, 1))
nc_close(newET24hNCDF4)

print(proc.time())
