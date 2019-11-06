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

nc<-nc_open(paste(fic.preproc, "/", fic, "_alb.nc", sep=""), write=TRUE, readunlim=FALSE, verbose=TRUE, auto_GMT=FALSE, suppress_dimvals=FALSE)

# Getting lat and lon values from old NetCDF
oldLat <- ncvar_get(nc, "lat", start=1, count=raster.elevation@nrows)
oldLon <- ncvar_get(nc, "lon", start=1, count=raster.elevation@ncols)

# Defining latitude and longitude dimensions
dimLatDef <- ncdim_def("lat", "degrees", oldLat, unlim=FALSE, longname="latitude")
dimLonDef <- ncdim_def("lon", "degrees", oldLon, unlim=FALSE, longname="longitude")

nc_close(nc)

alb <- raster(paste(fic.preproc, "/", fic, "_alb.nc", sep=""))
TS <- raster(paste(fic.preproc, "/", fic, "_TS.nc", sep=""))
NDVI <- raster(paste(fic.preproc, "/", fic, "_NDVI.nc", sep=""))
G <- raster(paste(fic.preproc, "/", fic, "_G.nc", sep=""))
Rn <- raster(paste(fic.preproc, "/", fic, "_Rn.nc", sep=""))

#################################### Constants ##########################################

k <- 0.41		# Von K?rm?n
g <- 9.81		# Gravity
rho <- 1.15		# Air density
cp <- 1004		# Specific heat of air
Gsc <- 0.082		# Solar constant (0.0820 MJ m-2 min-1)

proc.time()

phase2 <- function() {
  ######################## Selection of reference pixels ###################################
  
  # Hot Pixel Candidates
  HO<-Rn[]-G[] # Read as a Vector
  x<-TS[][(NDVI[]>0.15 &!is.na(NDVI[]))  & (NDVI[]<0.20 &!is.na(NDVI[])) ] # Returns a vector
  x<-x[x>273.16]
  TS.c.hot<-sort(x)[round(0.99*length(x))] # Returns one value
  HO.c.hot<-HO[][(NDVI[]>0.15 &!is.na(NDVI[])) & (NDVI[]<0.20 &!is.na(NDVI[])) & TS[]==TS.c.hot] # Returns one value
  
  print (length(HO.c.hot))

  if (length(HO.c.hot)==1){
    ll.hot<-which(TS[]==TS.c.hot & HO[]==HO.c.hot)
    xy.hot <- xyFromCell(TS, ll.hot)
    ll.hot.f<-cbind(as.vector(xy.hot[1,1]), as.vector(xy.hot[1,2]))
  }else{
    HO.c.hot.min<-sort(HO.c.hot)[ceiling(0.25*length(HO.c.hot))]
    HO.c.hot.max<-sort(HO.c.hot)[round(0.75*length(HO.c.hot))]
    ll.hot<-which(TS[]==TS.c.hot & HO[]>=HO.c.hot.min & HO[]<=HO.c.hot.max)
    xy.hot <- xyFromCell(TS, ll.hot)
    NDVI.hot<-extract(NDVI,xy.hot, buffer=105)
    NDVI.hot.2<-NDVI.hot[!sapply(NDVI.hot, is.null)]
    NDVI.hot.cv <- sapply(NDVI.hot.2,sd, na.rm=TRUE)/sapply(NDVI.hot.2, mean, na.rm=TRUE)
    i.NDVI.hot.cv<-which.min(NDVI.hot.cv)
    ll.hot.f<-cbind(as.vector(xy.hot[i.NDVI.hot.cv,1]), as.vector(xy.hot[i.NDVI.hot.cv,2]))
  }
  
  # Cold Pixel Candidates
  x<-TS[][(NDVI[]<0 &!is.na(NDVI[])) & !is.na(HO)]
  x<-x[x>273.16]
  TS.c.cold<-sort(x)[round(0.5*length(x))]
  HO.c.cold<-HO[(NDVI[]<0 & !is.na(NDVI[])) & TS[]==TS.c.cold & !is.na(HO)]
  if (length(HO.c.cold)==1){
    ll.cold<-which(TS[]==TS.c.cold & HO==HO.c.cold)
    xy.cold <- xyFromCell(TS, ll.cold)
    ll.cold.f<-cbind(as.vector(xy.cold[1,1]), as.vector(xy.cold[1,2]))
  }else{
    HO.c.cold.min<-sort(HO.c.cold)[ceiling(0.25*length(HO.c.cold))]
    HO.c.cold.max<-sort(HO.c.cold)[round(0.75*length(HO.c.cold))]
    ll.cold<-which(TS[]==TS.c.cold & (HO>=HO.c.cold.min &!is.na(HO)) & (HO<=HO.c.cold.max & !is.na(HO)))
    xy.cold <- xyFromCell(TS, ll.cold)
    NDVI.cold<-extract(NDVI,xy.cold, buffer=105)
    NDVI.cold.2<-NDVI.cold[!sapply(NDVI.cold, is.null)]
    
    # Maximum number of neighboring pixels with $NVDI < 0$
    t<-function(x){ sum(x<0,na.rm = TRUE)}
    n.neg.NDVI<-sapply(NDVI.cold.2,t)
    i.NDVI.cold<-which.max(n.neg.NDVI)
    ll.cold.f<-cbind(as.vector(xy.cold[i.NDVI.cold,1]), as.vector(xy.cold[i.NDVI.cold,2]))
  }
  
  # Location of reference pixels (hot and cold)
  ll_ref<-rbind(ll.hot.f[1,],ll.cold.f[1,])
  colnames(ll_ref)<-c("long", "lat")
  rownames(ll_ref)<-c("hot","cold")
  
  print(proc.time())
  
  ####################################################################################
  
  # Weather station data
  x<-3 					# Wind speed sensor Height (meters)
  hc<-0.2 				# Vegetation height (meters)
  Lat<-table.sw$V4[1] 	# Station Latitude
  Long<-table.sw$V5[1] 	# Station Longitude
  
  # Surface roughness parameters in station
  zom.est <- hc*0.12
  azom <- -3		#Parameter for the Zom image
  bzom <- 6.47	#Parameter for the Zom image
  F_int <- 0.16	#internalization factor for Rs 24 calculation (default value)
  
  print(proc.time())
  
  # Friction velocity at the station (ustar.est)
  ustar.est<-k*table.sw$V6[hour.image.station]/log((x)/zom.est)#Ajustar
  
  # Velocity 200 meters
  u200<-ustar.est/k*log(200/zom.est)
  
  # Zom for all pixels
  zom<-exp(azom+bzom*NDVI[]) # Changed from Raster to Vector
  
  # Initial values
  ustar<-NDVI
  ustar[]<-k*u200/(log(200/zom))			# Friction velocity for all pixels #RASTER - VETOR 
  rah<-NDVI
  rah[]<-(log(2/0.1))/(ustar[]*k) 		# Aerodynamic resistance for all pixels #RASTER - VETOR
  base_ref<-stack(NDVI,TS,Rn,G,ustar,rah) # Raster
  nbase<-c("NDVI","TS","Rn","G")
  names(base_ref)<-c(nbase,"ustar","rah")
  value.pixels.ref<-extract(base_ref,ll_ref)
  rownames(value.pixels.ref)<-c("hot","cold")
  H.hot<-value.pixels.ref["hot","Rn"]-value.pixels.ref["hot","G"]  
  value.pixel.rah<-value.pixels.ref["hot","rah"]
  
  i<-1
  Erro<-TRUE
  
  print(proc.time())
  
  # Beginning of the cycle stability
  while(Erro){
    rah.hot.0<-value.pixel.rah[i] # Value
    
    # Hot and cold pixels      
    dt.hot<-H.hot*rah.hot.0/(rho*cp) # Value
    b<-dt.hot/(value.pixels.ref["hot","TS"]-value.pixels.ref["cold","TS"]) # Value
    a<- -b*(value.pixels.ref["cold","TS"]-273.15) # Value
    
    # All pixels
    H<-rho*cp*(a+b*(TS[]-273.15))/rah[] # Changed from Raster to Vector
    L<- -1*((rho*cp*ustar[]^3*TS[])/(k*g*H)) # Changed from Raster to Vector
    y_0.1<-(1-16*0.1/L)^0.25 # Changed from Raster to Vector
    y_2<-(1-16*2/L)^0.25 # Changed from Raster to Vector
    x200<-(1-16*200/L)^0.25 # Changed from Raster to Vector
    psi_0.1<-2*log((1+y_0.1^2)/2) # Changed from Raster to Vector
    psi_0.1[L>0 &!is.na(L)]<--5*(0.1/L[L>0 &!is.na(L)]) # Changed from Raster to Vector
    psi_2<-2*log((1+y_2^2)/2)  # Changed from Raster to Vector
    psi_2[L>0 &!is.na(L) ]<--5*(2/L[L>0 &!is.na(L)]) # Changed from Raster to Vector
    psi_200<-2*log((1+x200)/2)+log((1+x200^2)/2)-2*atan(x200)+0.5*pi # Changed from Raster to Vector
    psi_200[L>0 &!is.na(L) ]<--5*(2/L[(L>0 &!is.na(L))]) # Changed from Raster to Vector
    ustar<-k*u200/(log(200/zom)-psi_200) # Changed from Raster to Vector # Friction velocity for all pixels
    rah<-NDVI
    rah[]<-(log(2/0.1)-psi_2+psi_0.1)/(ustar*k) # Changed from Raster to Vector # Aerodynamic resistency for all pixels
    rah.hot<-extract(rah,matrix(ll_ref["hot",],1,2)) # Value
    value.pixel.rah<-c(value.pixel.rah,rah.hot) # Value
    Erro<-(abs(1-rah.hot.0/rah.hot)>=0.1)
    i<-i+1
  }
  
  print(proc.time())
  
  # End sensible heat flux (H)
  
  # Hot and cold pixels
  dt.hot<-H.hot*rah.hot/(rho*cp)                  
  b<-dt.hot/(value.pixels.ref["hot","TS"]-value.pixels.ref["cold","TS"]) 
  a<- -b*(value.pixels.ref["cold","TS"]-273.15)                          
  
  print(proc.time())
  
  # All pixels
  
  H<-rho*cp*(a+b*(TS[]-273.15))/rah[] # Vector
  subset <- !is.na(H) & H > (Rn[]-G[])
  H[subset]<-(Rn[]-G[])*which(subset) # Vector
  
  print(proc.time())
  
  # Instant latent heat flux (LE)
  LE<-Rn[]-G[]-H
  
  # Upscalling temporal
  dr<-(1/dist.dia.juliano)^2 		# Inverse square of the distance on Earth-SOL
  sigma<-0.409*sin(((2*pi/365)*Dia.juliano)-1.39) # Declination Solar (rad)
  phi<-(pi/180)*Lat 								# Solar latitude in degrees
  omegas<-acos(-tan(phi)*tan(sigma)) 				# Angle Time for sunsets (rad)
  Ra24h<-(((24*60/pi)*Gsc*dr)*(omegas*sin(phi)*
                                 sin(sigma)+cos(phi)*cos(sigma)*sin(omegas)))*(1000000/86400)
  
  print(proc.time())
  
  # Short wave radiation incident in 24 hours (Rs24h)
  Rs24h<-F_int*sqrt(max(table.sw$V7[])-min(table.sw$V7[]))*Ra24h
  
  FL<-110                                
  Rn24h_dB<-(1-alb[])*Rs24h-FL*Rs24h/Ra24h		# Method of Bruin #VETOR
  
  # Evapotranspiration fraction Bastiaanssen
  EF<-NDVI
  EF[]<-LE/(Rn[]-G[])
  
  # Sensible heat flux 24 hours (H24h)
  H24h_dB<-(1-EF[])*Rn24h_dB
  
  # Latent Heat Flux 24 hours (LE24h)
  LE24h_dB<-EF[]*Rn24h_dB
  
  # Evapotranspiration 24 hours (ET24h)
  ET24h_dB<-NDVI
  ET24h_dB[]<-LE24h_dB*86400/((2.501-0.00236* (max(table.sw$V7[])+min(table.sw$V7[]))/2)*10^6)
  
  print(proc.time())
  
  output.evapo<-stack(EF,ET24h_dB)
  names(output.evapo)<-c('EF','ET24h')
  writeRaster(output.evapo, output.path, overwrite=TRUE, format="CDF", varname=fic, varunit="daily", longname=fic, xname="lon", yname="lat", bylayer=TRUE, suffix="names")
  
  print(proc.time())
  
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
}

tryCatch({
  res <- withTimeout({
    phase2();
  }, timeout=5400);
}, TimeoutException=function(ex) {
  cat("Image phase two processing timedout. Exiting with 124 code...\n");
  quit("no", 124, FALSE)
})

proc.time()

