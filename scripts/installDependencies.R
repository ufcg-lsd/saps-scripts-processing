if (!require(devtools)) {
  install.packages("devtools", repos = "http://cran.us.r-project.org")
}

if (packageVersion("devtools") < 1.6) {
  install.packages("devtools", repos = "http://cran.us.r-project.org")
}

if(!require(sp)){
  install.packages("sp", repos = "http://cran.us.r-project.org")
}

if(!require(rgdal)){
  install.packages("rgdal", repos = "http://cran.us.r-project.org")
}

if(!require(raster)){
  install.packages("raster", repos = "http://cran.us.r-project.org")
}

if(!require(gstat)){
  install.packages("gstat", repos = "http://cran.us.r-project.org")
}

if(!require(httr)){
  install.packages("httr", repos = "http://cran.us.r-project.org")
}

# Install in shell R using command "devtools::install_github("gowusu/sebkc", force = TRUE)"

#if(!require(sebkc)){
#   devtools::install_github("gowusu/sebkc", force = TRUE)
#}

if (!require(R.utils)) {
  install.packages("R.utils", repos = "http://cran.us.r-project.org")
}

if (!require(maptools)) {
  install.packages("maptools", repos = "http://cran.us.r-project.org")
}

if (!require(ncdf4)) {
  install.packages("ncdf4", repos = "http://cran.us.r-project.org")
}

if (!require(snow)) {
  install.packages("snow", repos = "http://cran.us.r-project.org")
}
