#' Spatial prediction data setup
#' M. Walsh, December 2016

# install.packages(c("downloader","rgdal","raster"), dependencies=T)
require(downloader)
require(rgdal)
require(raster)

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("Spatial_data", showWarnings=F)
setwd("./Spatial_data")

# Download GeoSurvey data
download("https://www.dropbox.com/s/6gr2s8s7zgj7hr7/GS_CNG_2016.csv.zip?dl=0", "GS_CNG_2016.csv.zip", mode="wb")
unzip("GS_CNG_2016.csv.zip", overwrite=T)
geo <- read.table("GS_CNG_2016.csv", header=T, sep=",")

# Download MobileSurvey data
download("https://www.dropbox.com/s/lbzcnzrqbuiyh8r/MS_CNG_2016.csv.zip?dl=0", "MS_CNG_2016.csv.zip", mode="wb")
unzip("MS_CNG_2016.csv.zip", overwrite=T)
mob <- read.table("MS_CNG_2016.csv", header=T, sep=",")

# Download MIR prediction data
download("https://www.dropbox.com/s/kjm69pdsdc21hxg/Top_MIR_pred.csv.zip?dl=0", "Top_MIR_pred.csv.zip", mode="wb")
unzip("Top_MIR_pred.csv.zip", overwrite=T)
mir <- read.table("Top_MIR_pred.csv", header=T, sep=",")

# Download grids
download("https://www.dropbox.com/s/vewvi0l1o949yh2/OCP_grids.zip?dl=0", "OCP_grids.zip", mode="wb")
unzip("OCP_grids.zip", overwrite=T)
glist <- list.files(pattern="tif", full.names=T)
grids <- stack(glist)

# Overlay surveys with gridded covariates ---------------------------------
# Extract grids at GeoSurvey locations
geo.proj <- as.data.frame(project(cbind(geo$Lon, geo$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(geo.proj) <- c("x","y") ## laea coordinates
geo <- cbind(geo, geo.proj)
coordinates(geo) <- ~x+y
projection(geo) <- projection(grids)
geogrd <- extract(grids, geo)
geogrd <- as.data.frame(geogrd)
geo <- cbind.data.frame(geo, geogrd)
geo <- unique(na.omit(geo)) ## includes only unique & complete records

# Extract grids at MobileSurvey locations
mob.proj <- as.data.frame(project(cbind(mob$Lon, mob$Lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(mob.proj) <- c("x","y") ## laea coordinates
mob <- cbind(mob, mob.proj)
coordinates(mob) <- ~x+y
projection(mob) <- projection(grids)
mobgrd <- extract(grids, mob)
mobgrd <- as.data.frame(mobgrd)
mob <- cbind.data.frame(mob, mobgrd)
mob <- unique(na.omit(mob)) ## includes only unique & complete records

# Extract grids at MIR locations
coordinates(mir) <- ~x+y
projection(mir) <- projection(grids)
mirgrd <- extract(grids, mir)
mirgrd <- as.data.frame(mirgrd)
mir <- cbind(mir, mirgrd)
mir <- unique(na.omit(mir)) ## includes only unique & complete records

# Train/Test set partitions -----------------------------------------------




