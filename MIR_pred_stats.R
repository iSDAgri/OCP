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
C <-  read.table("C_pred.csv", header=T, sep=",")
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

# Sample ID codebook
iita <- read.table("IITA_codebook.csv", header=T, sep=",")
ssid <- read.table("OCP_SSID.csv", header=T, sep=",")
ssid <- merge(ssid, iita, by="ssid")

# AfSIS reference data
download("https://www.dropbox.com/s/srhas5aa6b76kb3/AfSIS_ref_data.csv.zip?dl=0", "AfSIS_ref_data.csv.zip", mode="wb")
unzip("AfSIS_ref_data.csv.zip", overwrite=T)
ref <- read.table("AfSIS_ref_data.csv", header=T, sep=",")
ref <- ref[which(ref$Depth==10), ] ## select topsoils
ref <- na.omit(ref) ## omit any missing values

# Topsoil prediction data setup -------------------------------------------
# Topsoil Nitrogen predictions
top_N <- merge(ssid, N, by="SSN")
top_N <- top_N[which(top_N$depth=="top"), ]
top_N <- top_N[!duplicated(top_N[,2]), ]

# Topsoil Phosporus predictions
top_P <- merge(ssid, P, by="SSN")
top_P <- top_P[which(top_P$depth=="top"), ]
top_P <- top_P[!duplicated(top_P[,2]), ]

# Topsoil Potassium predictions
top_K <- merge(ssid, K, by="SSN")
top_K <- top_K[which(top_K$depth=="top"), ]
top_K <- top_K[!duplicated(top_K[,2]), ]

# Topsoil Sulfur predictions
top_S <- merge(ssid, S, by="SSN")
top_S <- top_S[which(top_S$depth=="top"), ]
top_S <- top_S[!duplicated(top_S[,2]), ]

# Topsoil Boron predictions
top_B <- merge(ssid, B, by="SSN")
top_B <- top_B[which(top_B$depth=="top"), ]
top_B <- top_B[!duplicated(top_B[,2]), ]

# Topsoil Zinc predictions
top_Zn <- merge(ssid, Zn, by="SSN")
top_Zn <- top_Zn[which(top_Zn$depth=="top"), ]
top_Zn <- top_Zn[!duplicated(top_Zn[,2]), ]

# Topsoil cum distribution plots ------------------------------------------
par(mfrow=c(3,2), mar=c(5,4.5,1,1))
plo <- 0.25
phi <- 0.75

# Topsoil Nitrogen predictions
plot(ecdf(top_N$p50), main="", xlab="N (%)", cex.lab=1.4, ylab="CDF", xlim=c(0,0.3), verticals=T, lty=1, lwd=2, col="red", do.points=F, col.01line = NULL)
abline(v=quantile(ref$N/10000, probs=plo), h=plo, lty=1, col="grey")
abline(v=quantile(ref$N/10000, probs=phi), h=phi, lty=1, col="grey")
plot(ecdf(top_N$p05), add=T, verticals=T, lty=2, lwd=2, col="dark grey", do.points=F, col.01line = NULL)
plot(ecdf(top_N$p95), add=T, verticals=T, lty=2, lwd=2, col="dark grey", do.points=F, col.01line = NULL)

# Topsoil Phosphorus predictions
plot(ecdf(top_P$p50), main="", xlab="P (ppm)", cex.lab=1.4, ylab="CDF", xlim=c(0,40), verticals=T, lty=1, lwd=2, col="red", do.points=F, col.01line = NULL)
abline(v=quantile(ref$P, probs=plo), h=plo, lty=1, col="grey")
abline(v=quantile(ref$P, probs=phi), h=phi, lty=1, col="grey")
plot(ecdf(top_P$p05), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)
plot(ecdf(top_P$p95), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)

# Topsoil Potassium predictions
plot(ecdf(top_K$p50), main="", xlab="K (ppm)", cex.lab=1.4, ylab="CDF", xlim=c(0,300), verticals=T, lty=1, lwd=2, col="red", do.points=F, col.01line = NULL)
abline(v=quantile(ref$K, probs=plo), h=plo, lty=1, col="grey")
abline(v=quantile(ref$K, probs=phi), h=phi, lty=1, col="grey")
plot(ecdf(top_K$p05), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)
plot(ecdf(top_K$p95), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)

# Topsoil Sulfur predictions
plot(ecdf(top_S$p50), main="", xlab="S (ppm)", cex.lab=1.4, ylab="CDF", xlim=c(0,40), verticals=T, lty=1, lwd=2, col="red", do.points=F, col.01line = NULL)
abline(v=quantile(ref$S, probs=plo), h=plo, lty=1, col="grey")
abline(v=quantile(ref$S, probs=phi), h=phi, lty=1, col="grey")
plot(ecdf(top_S$p05), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)
plot(ecdf(top_S$p95), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)

# Topsoil Boron predictions
plot(ecdf(top_B$p50), main="", xlab="B (ppm)", cex.lab=1.4, ylab="CDF", xlim=c(0,1), verticals=T, lty=1, lwd=2, col="red", do.points=F, col.01line = NULL)
abline(v=quantile(ref$B, probs=plo), h=plo, lty=1, col="grey")
abline(v=quantile(ref$B, probs=phi), h=phi, lty=1, col="grey")
plot(ecdf(top_B$p05), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)
plot(ecdf(top_B$p95), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)

# Topsoil Zinc predictions
plot(ecdf(top_Zn$p50), main="", xlab="Zn (ppm)", cex.lab=1.4, ylab="CDF", xlim=c(0,4), verticals=T, lty=1, lwd=2, col="red", do.points=F, col.01line = NULL)
abline(v=quantile(ref$Zn, probs=plo), h=plo, lty=1, col="grey")
abline(v=quantile(ref$Zn, probs=phi), h=phi, lty=1, col="grey")
plot(ecdf(top_Zn$p05), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)
plot(ecdf(top_Zn$p95), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)

par(mfrow=c(1,1))

# MIR prediction data reshape ---------------------------------------------
# Topsoil Nitrogen predictions
names(top_N)[58:60] <- c("p.5", "p.50", "p.95")
top_N <- reshape(top_N, direction="long", varying=58:60, idvar="ssid", v.names="N", timevar="plevel") ## long format
top_N$lo <- ifelse(top_N$N < quantile(ref$N/10000, probs=plo), 1, 0) ## identifies low levels
top_N$hi <- ifelse(top_N$N > quantile(ref$N/10000, probs=phi), 1, 0) ## identifies high levels

# Topsoil Phosphorus predictions
names(top_P)[58:60] <- c("p.5", "p.50", "p.95")
top_P <- reshape(top_P, direction="long", varying=58:60, idvar="ssid", v.names="P", timevar="plevel") ## long format
top_P$lo <- ifelse(top_P$P < quantile(ref$P, probs=plo), 1, 0) ## identifies low levels
top_P$hi <- ifelse(top_P$P > quantile(ref$P, probs=phi), 1, 0) ## identifies high levels

# Topsoil Potassium predictions
names(top_K)[58:60] <- c("p.5", "p.50", "p.95")
top_K <- reshape(top_K, direction="long", varying=58:60, idvar="ssid", v.names="P", timevar="plevel") ## long format
top_K$lo <- ifelse(top_K$K < quantile(ref$K, probs=plo), 1, 0) ## identifies low levels
top_K$hi <- ifelse(top_K$K > quantile(ref$K, probs=phi), 1, 0) ## identifies high levels

# Topsoil Sulfur predictions
names(top_S)[58:60] <- c("p.5", "p.50", "p.95")
top_S <- reshape(top_S, direction="long", varying=58:60, idvar="ssid", v.names="P", timevar="plevel") ## long format
top_S$lo <- ifelse(top_S$S < quantile(ref$S, probs=plo), 1, 0) ## identifies low levels
top_S$hi <- ifelse(top_S$S > quantile(ref$S, probs=phi), 1, 0) ## identifies high levels

# Topsoil Boron predictions
names(top_B)[58:60] <- c("p.5", "p.50", "p.95")
top_B <- reshape(top_B, direction="long", varying=58:60, idvar="ssid", v.names="P", timevar="plevel") ## long format
top_B$lo <- ifelse(top_B$B < quantile(ref$B, probs=plo), 1, 0) ## identifies low levels
top_B$hi <- ifelse(top_B$B > quantile(ref$B, probs=phi), 1, 0) ## identifies high levels

# Topsoil Zinc predictions
names(top_Zn)[58:60] <- c("p.5", "p.50", "p.95")
top_Zn <- reshape(top_Zn, direction="long", varying=58:60, idvar="ssid", v.names="P", timevar="plevel") ## long format
top_Zn$lo <- ifelse(top_Zn$Zn < quantile(ref$Zn, probs=plo), 1, 0) ## identifies low levels
top_Zn$hi <- ifelse(top_Zn$Zn > quantile(ref$Zn, probs=phi), 1, 0) ## identifies high levels
