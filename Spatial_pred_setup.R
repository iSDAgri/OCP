#' Central Nigeria GeoSurvey, MobileSurvey & soil data spatial prediction setup 
#' M. Walsh, December 2016

# install.packages(c("downloader","rgdal","raster"), dependencies=T)
suppressPackageStartupMessages({
  require(downloader)
  require(rgdal)
  require(raster)
})

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("GRID_data", showWarnings=F)
setwd("./GRID_data")

# Download GeoSurvey data
download("https://www.dropbox.com/s/6gr2s8s7zgj7hr7/GS_CNG_2016.csv.zip?raw=1", "GS_CNG_2016.csv.zip", mode="wb")
unzip("GS_CNG_2016.csv.zip", overwrite=T)
geo <- read.table("GS_CNG_2016.csv", header=T, sep=",")

# Download MobileSurvey data
download("https://www.dropbox.com/s/lbzcnzrqbuiyh8r/MS_CNG_2016.csv.zip?raw=1", "MS_CNG_2016.csv.zip", mode="wb")
unzip("MS_CNG_2016.csv.zip", overwrite=T)
mob <- read.table("MS_CNG_2016.csv", header=T, sep=",")

# Download MIR prediction data
download("https://www.dropbox.com/s/x2or21x624n75gu/Top_MIR_pred.csv.zip?raw=1", "Top_MIR_pred.csv.zip", mode="wb")
unzip("Top_MIR_pred.csv.zip", overwrite=T)
mir <- read.table("Top_MIR_pred.csv", header=T, sep=",")

# Download grids
download("https://www.dropbox.com/s/vewvi0l1o949yh2/OCP_grids.zip?raw=1", "OCP_grids.zip", mode="wb")
unzip("OCP_grids.zip", overwrite=T)
download("https://www.dropbox.com/s/03sl7rnfpdb2d1o/LGA_ID.zip?raw=1", "LGA_ID.zip", mode="wb")
unzip("LGA_ID.zip", overwrite=T)
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
mirvars <- c("ssid","x","y","N","P","K","S","B","Zn") ## select variables
nut <- mir[mirvars]
coordinates(nut) <- ~x+y
projection(nut) <- projection(grids)
nutgrd <- extract(grids, nut)
nutgrd <- as.data.frame(nutgrd)
nut <- cbind.data.frame(nut, nutgrd)
nut <- unique(na.omit(nut)) ## includes only unique & complete records

# extract variables for nutrient requirement summaries
reqvars <- c("ssid","x","y","LGA_ID","N","P","K","S","B","Zn","BD20")
req <- nut[reqvars]

# Write files -------------------------------------------------------------
write.csv(geo, "geo.csv", row.names = F)
write.csv(mob, "mob.csv", row.names = F)
write.csv(nut, "nut.csv", row.names = F)
write.csv(req, "req.csv", row.names = F)

# remove extraneous objects from memory
rm(list=setdiff(ls(), c("geo", "mob", "nut", "req", "grids")))
