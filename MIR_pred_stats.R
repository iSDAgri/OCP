#' Data summaries of MIR-M3 prediction results
#' M. Walsh, June 2016

# install.packages(c("downloader","compositions"), dependencies=T)
require(downloader)

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("MIR_data", showWarnings=F)
setwd("./MIR_data")

# Download M3 prediction data
download("https://www.dropbox.com/s/uyelu38ah6zos6h/OCP_KBR_pred.zip?dl=0", "OCP_KBR_pred.zip", mode="wb")
unzip("OCP_KBR_pred.zip", overwrite=T)
N <-  read.table("N_pred.csv", header=T, sep=",")
B <-  read.table("B_pred.csv", header=T, sep=",")
Mg <- read.table("Mg_pred.csv", header=T, sep=",")
P <-  read.table("P_pred.csv", header=T, sep=",")
S <-  read.table("S_pred.csv", header=T, sep=",")
K <-  read.table("K_pred.csv", header=T, sep=",")
Ca <- read.table("Ca_pred.csv", header=T, sep=",")
Mn <- read.table("Mn_pred.csv", header=T, sep=",")
Fe <- read.table("Fe_pred.csv", header=T, sep=",")
Cu <- read.table("Cu_pred.csv", header=T, sep=",")
Zn <- read.table("Zn_pred.csv", header=T, sep=",")

# AfSIS reference data
download("https://www.dropbox.com/s/srhas5aa6b76kb3/AfSIS_ref_data.csv.zip?dl=0", "AfSIS_ref_data.csv.zip", mode="wb")
unzip("AfSIS_ref_data.csv.zip", overwrite=T)
ref <- read.table("AfSIS_ref_data.csv", header=T, sep=",")
ref$N <- ref$N/10000 ## convert from ppm to %
ref <- ref[which(ref$Depth==10), ] ## select topsoils
ref <- na.omit(ref) ## omit any missing values

# Cum distribution plots --------------------------------------------------
par(mfrow=c(3,2), mar=c(5,4.5,1,1))
plo <- 0.25
phi <- 0.75

# Nitrogen
plot(ecdf(N$p50), main="", xlab="N (%)", cex.lab=1.4, ylab="CDF", xlim=c(0,0.3), verticals=T, lty=1, lwd=2, col="red", do.points=F, col.01line = NULL)
abline(v=quantile(ref$N, probs=plo), h=plo, lty=1, col="grey")
abline(v=quantile(ref$N, probs=phi), h=phi, lty=1, col="grey")
plot(ecdf(N$p05), add=T, verticals=T, lty=2, lwd=2, col="dark grey", do.points=F, col.01line = NULL)
plot(ecdf(N$p95), add=T, verticals=T, lty=2, lwd=2, col="dark grey", do.points=F, col.01line = NULL)

# Phosphorus
plot(ecdf(P$p50), main="", xlab="P (ppm)", cex.lab=1.4, ylab="CDF", xlim=c(0,40), verticals=T, lty=1, lwd=2, col="red", do.points=F, col.01line = NULL)
abline(v=quantile(ref$P, probs=plo), h=plo, lty=1, col="grey")
abline(v=quantile(ref$P, probs=phi), h=phi, lty=1, col="grey")
plot(ecdf(P$p05), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)
plot(ecdf(P$p95), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)

# Potassium
plot(ecdf(K$p50), main="", xlab="K (ppm)", cex.lab=1.4, ylab="CDF", xlim=c(0,300), verticals=T, lty=1, lwd=2, col="red", do.points=F, col.01line = NULL)
abline(v=quantile(ref$K, probs=plo), h=plo, lty=1, col="grey")
abline(v=quantile(ref$K, probs=phi), h=phi, lty=1, col="grey")
plot(ecdf(K$p05), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)
plot(ecdf(K$p95), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)

# Sulfur
plot(ecdf(S$p50), main="", xlab="S (ppm)", cex.lab=1.4, ylab="CDF", xlim=c(0,40), verticals=T, lty=1, lwd=2, col="red", do.points=F, col.01line = NULL)
abline(v=quantile(ref$S, probs=plo), h=plo, lty=1, col="grey")
abline(v=quantile(ref$S, probs=phi), h=phi, lty=1, col="grey")
plot(ecdf(S$p05), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)
plot(ecdf(S$p95), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)

# Boron
plot(ecdf(B$p50), main="", xlab="B (ppm)", cex.lab=1.4, ylab="CDF", xlim=c(0,1), verticals=T, lty=1, lwd=2, col="red", do.points=F, col.01line = NULL)
abline(v=quantile(ref$B, probs=plo), h=plo, lty=1, col="grey")
abline(v=quantile(ref$B, probs=phi), h=phi, lty=1, col="grey")
plot(ecdf(B$p05), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)
plot(ecdf(B$p95), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)

# Zinc
plot(ecdf(Zn$p50), main="", xlab="Zn (ppm)", cex.lab=1.4, ylab="CDF", xlim=c(0,4), verticals=T, lty=1, lwd=2, col="red", do.points=F, col.01line = NULL)
abline(v=quantile(ref$Zn, probs=plo), h=plo, lty=1, col="grey")
abline(v=quantile(ref$Zn, probs=phi), h=phi, lty=1, col="grey")
plot(ecdf(Zn$p05), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)
plot(ecdf(Zn$p95), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)
