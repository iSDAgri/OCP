#' Mehlich-3 (M3) data checks by comparison to CNLS reanalyses
#' M. Walsh, March 2017

# install.packages(c("downloader","compositions"), dependencies=T)
suppressPackageStartupMessages({
  require(downloader)
})

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("M3_data", showWarnings=F)
setwd("./M3_data")
getwd() ## your current working directory

# Data download
download("https://www.dropbox.com/s/kh8vewc9kas3x8k/OCP_M3_data.zip?raw=1", "OCP_M3_data.zip", mode="wb")
unzip("OCP_M3_data.zip", overwrite=T)
cbook <- read.table("SSID.csv", header=T, sep=",")  ## read lab soil ID codebook
iita1 <- read.table("iita1.csv", header=T, sep=",") ## IITA M3 initial data in mg/kg
iita2 <- read.table("iita2.csv", header=T, sep=",") ## IITA M3 rerun data in mg/kg
cnlsr <- read.table("cnls.csv", header=T, sep=",")  ## CNLS ICP-MS M3 reference data
cnlsr <- merge(cbook, cnlsr, by="SSN")

# Assemble dataframes
cnlsr$Ca <- cnlsr$ExCa * 200 ## ExCa conversion to ppm
cnlsr$Mg <- cnlsr$ExMg * 120 ## ExMg conversion to ppm
cnlsr$K  <- cnlsr$ExK * 391  ## ExK conversion to ppm

iita2 <- iita2[c(1,5,7,6,8,4,3,9,11,12,10)]
cnlsr <- cnlsr[c(3,10,23,11,21:22,6,9,7,12,8)]
names(cnlsr) <- c("SSID","P","K","S","Ca","Mg","B","Mn","Cu","Zn","Fe")
compr <- merge(cnlsr, iita2, by="SSID")

# Plots -------------------------------------------------------------------
# Macro-nutrients
par(mfrow=c(2,2), mar=c(5,5,1,1))
plot(P.y~P.x, cex=1.2, xlab=expression(paste("P"[r], " (ppm)")), ylab=expression(paste("P"[s], " (ppm)")), 
     cex.lab=1.5, xlim=c(0, max(compr$P.x)), ylim= c(0, max(compr$P.x)), compr)
abline(0,1)
plot(K.y~K.x, cex=1.2, xlab=expression(paste("K"[r], " (ppm)")), ylab=expression(paste("K"[s], " (ppm)")), 
     cex.lab=1.5, xlim=c(0, max(compr$K.x)), ylim= c(0, max(compr$K.x)), compr)
abline(0,1)
plot(S.y~S.x, cex=1.2, xlab=expression(paste("S"[r], " (ppm)")), ylab=expression(paste("S"[s], " (ppm)")), 
     cex.lab=1.5, xlim=c(0, max(compr$S.x)), ylim= c(0, max(compr$S.x)), compr)
abline(0,1)
plot(Ca.y~Ca.x, cex=1.2, xlab=expression(paste("Ca"[r], " (ppm)")), ylab=expression(paste("Ca"[s], " (ppm)")), 
     cex.lab=1.5, xlim=c(0, max(compr$Ca.x)), ylim= c(0, max(compr$Ca.x)), compr)
abline(0,1)

# Micro-nutrients
par(mfrow=c(2,2), mar=c(5,5,1,1))
plot(B.y~B.x, cex=1.2, xlab=expression(paste("B"[r], " (ppm)")), ylab=expression(paste("B"[s], " (ppm)")), 
     cex.lab=1.5, xlim=c(0, max(compr$B.x)), ylim= c(0, max(compr$B.x)), compr)
abline(0,1)
plot(Mn.y~Mn.x, cex=1.2, xlab=expression(paste("Mn"[r], " (ppm)")), ylab=expression(paste("Mn"[s], " (ppm)")), 
     cex.lab=1.5, xlim=c(0, max(compr$Mn.x)), ylim= c(0, max(compr$Mn.x)), compr)
abline(0,1)
plot(Cu.y~Cu.x, cex=1.2, xlab=expression(paste("Cu"[r], " (ppm)")), ylab=expression(paste("Cu"[s], " (ppm)")), 
     cex.lab=1.5, xlim=c(0, max(compr$Cu.x)), ylim= c(0, max(compr$Cu.x)), compr)
abline(0,1)
plot(Zn.y~Zn.x, cex=1.2, xlab=expression(paste("Zn"[r], " (ppm)")), ylab=expression(paste("Zn"[s], " (ppm)")), 
     cex.lab=1.5, xlim=c(0, max(compr$Zn.x)), ylim= c(0, max(compr$Zn.x)), compr)
abline(0,1)

# Mg & Fe
par(mfrow=c(1,2), mar=c(5,5,1,1))
plot(Mg.y~Mg.x, cex=1.2, xlab=expression(paste("Mg"[r], " (ppm)")), ylab=expression(paste("Mg"[s], " (ppm)")), 
     cex.lab=1.5, xlim=c(0, max(compr$Mg.x)), ylim= c(0, max(compr$Mg.x)), compr)
abline(0,1)
plot(Fe.y~Fe.x, cex=1.2, xlab=expression(paste("Fe"[r], " (ppm)")), ylab=expression(paste("Fe"[s], " (ppm)")), 
     cex.lab=1.5, xlim=c(0, max(compr$Fe.x)), ylim= c(0, max(compr$Fe.x)), compr)
abline(0,1)
