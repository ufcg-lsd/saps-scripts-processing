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
dist.dia.juliano <- as.numeric(MTL$V2[MTL$V1 == grep(pattern="EARTH_SUN_DISTANCE", MTL$V1, value=TRUE)])

output.path<-paste(dados$Path.Output[1], "/", fic, ".nc", sep="")

# Time image
acquired_date <- as.Date(MTL$V2[MTL$V1==grep(pattern="DATE_ACQUIRED", MTL$V1, value=TRUE)])
daysSince1970 <- as.numeric(acquired_date)
tdim <- ncdim_def("time", "days since 1970-1-1", daysSince1970, unlim=TRUE, create_dimvar=TRUE, "standard", "time")

fic.preproc <- dados$Path.Prepoc[1]

raster.elevation <- raster(paste(fic.preproc, "/elevation.tif", sep=""))
albedo <- raster(paste(fic.preproc, "/", fic, "_alb.nc", sep=""))
Ts <- raster(paste(fic.preproc, "/", fic, "_TS.nc", sep=""))
NDVI <- raster(paste(fic.preproc, "/", fic, "_NDVI.nc", sep=""))
LAI <- raster(paste(fic.preproc, "/", fic, "_LAI.nc", sep=""))
SAVI <- raster(paste(fic.preproc, "/", fic, "_SAVI.nc", sep=""))

##### Opening files ######

azom<- -3    #Parametro para o Zom imagem
bzom<- 6.47  #Parametro para o Zom imagem
zom<- exp(azom+bzom*NDVI)

# Upscalling temporal
Lat<-table.sw$V4[1] 	                            # Station Latitude
Gsc <- 0.082		                                # Solar constant (0.0820 MJ m-2 min-1)
F_int <- 0.16	                                    # Internalization factor for Rs 24 calculation (default value)
dr<-(1/dist.dia.juliano)^2 		                    # Inverse square of the distance on Earth-SOL
sigma<-0.409*sin(((2*pi/365)*Dia.juliano)-1.39)     # Declination Solar (rad)
phi<-(pi/180)*Lat 					                # Solar latitude in degrees
omegas<-acos(-tan(phi)*tan(sigma)) 			        # Angle Time for sunsets (rad)
Ra24h<-(((24*60/pi)*Gsc*dr)*(omegas*sin(phi)*
        sin(sigma)+cos(phi)*cos(sigma)*sin(omegas)))*(1000000/86400)

# Short wave radiation incident in 24 hours (Rs24h)
Rs24h<-F_int*sqrt(max(table.sw$V7[])-min(table.sw$V7[]))*Ra24h

FL<-110                                
Rn24h<-(1-albedo[])*Rs24h-FL*Rs24h/Ra24h

proc.time()

series=tseb(Ts=Ts, LAI=LAI, DOY=Dia.juliano, xyhot="full",
        xycold="full", albedo=albedo, Tmax=max(table.sw$V7[]),
        Tmin=min(table.sw$V7[]), RHmax=NULL, RHmin=NULL, zom=zom,
        NDVI=NDVI, SAVI=NULL, hc=3, DEM=raster.elevation, viewangle=0,
        sunelev=sun_elevation, welev=385, 
        u=table.sw$V6[hour.image.station], s=2, C=90, lapse=0.0065,
        Rn24 = Rn24h, zx=200, zomw=0.2, xPT=1.3, clump=1, fg=1, fc="auto",
        network="series", latitude=table.sw$V4[1], n=12, ETr = 10)

proc.time()

output.evapo<-stack(series$EF, series$ET24)
names(output.evapo)<-c('EF','ET24h')
writeRaster(output.evapo, output.path, overwrite=TRUE, format="CDF", varname=fic, varunit="daily", longname=fic, xname="lon", yname="lat", bylayer=TRUE, suffix="names")

proc.time()

# Opening old EF NetCDF
var_output<-paste(dados$Path.Output[1],"/",fic,"_EF.nc",sep="")
nc<-nc_open(var_output, write=TRUE,readunlim=FALSE,verbose=TRUE,auto_GMT=FALSE,suppress_dimvals=FALSE)

# Getting lat and lon values from old NetCDF
oldLat <- ncvar_get(nc, "lat", start=1, count=raster.elevation@nrows)
oldLon <- ncvar_get(nc, "lon", start=1, count=raster.elevation@ncols)

# Defining latitude and longitude dimensions
dimLatDef <- ncdim_def("lat", "degrees", oldLat, unlim=FALSE, longname="latitude")
dimLonDef <- ncdim_def("lon", "degrees", oldLon, unlim=FALSE, longname="longitude")

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
