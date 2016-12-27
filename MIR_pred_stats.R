#' Data summaries of MIR prediction results
#' M. Walsh, June 2016

# install.packages(c("downloader","rgdal"), dependencies=T)
require(downloader)
require(rgdal)

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("MIR_data", showWarnings=F)
setwd("./MIR_data")

# Download MIR prediction data
download("https://www.dropbox.com/s/uyelu38ah6zos6h/OCP_KBR_pred.zip?dl=0", "OCP_KBR_pred.zip", mode="wb")
unzip("OCP_KBR_pred.zip", overwrite=T)
pH <- read.table("pH_pred.csv", header=T, sep=",")
EC <- read.table("EC_pred.csv", header=T, sep=",")
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
Al <- read.table("Al_pred.csv", header=T, sep=",")
Na <- read.table("Na_pred.csv", header=T, sep=",")
Hp <- read.table("Hp_pred.csv", header=T, sep=",")

# IITA Sample ID codebook
iita <- read.table("IITA_codebook.csv", header=T, sep=",")
ssid <- read.table("OCP_SSID.csv", header=T, sep=",")
ssid <- merge(ssid, iita, by="ssid")

# AfSIS reference data
download("https://www.dropbox.com/s/7igjky8m1brphdp/AfSIS_ref_data.csv.zip?dl=0", "AfSIS_ref_data.csv.zip", mode="wb")
unzip("AfSIS_ref_data.csv.zip", overwrite=T)
ref <- read.table("AfSIS_ref_data.csv", header=T, sep=",")
ref <- ref[which(ref$Depth==10), ] ## select topsoils
ref <- na.omit(ref) ## omit any missing values

# Topsoil prediction data setup -------------------------------------------
# Topsoil pH
pH_lo <- 5.5 ## low reference level
pH_hi <- 8.5 ## high reference level
top_pH <- merge(ssid, pH, by="SSN") ## attaches codebook, lon/lat etc
top_pH <- top_pH[which(top_pH$depth=="top"), ] 
top_pH <- top_pH[!duplicated(top_pH[,2]), ]
top_pH$pH_lo <- rowSums(top_pH[,8:57]< pH_lo) ## calculates number of MCMC draws < pH_lo
top_pH$pH_hi <- rowSums(top_pH[,8:57]> pH_hi) ## calculates number of MCMC draws > pH_hi
colnames(top_pH)[59] <- "pH"

# Topsoil electrical conductivity (EC)
EC_lo <- quantile(ref$EC*1000, probs=0.25) ## low reference level
EC_hi <- quantile(ref$EC*1000, probs=0.75) ## high reference level
top_EC <- merge(ssid, EC, by="SSN")
top_EC <- top_EC[which(top_EC$depth=="top"), ]
top_EC <- top_EC[!duplicated(top_EC[,2]), ]
top_EC$EC_lo <- rowSums(top_EC[,8:57]< EC_lo)
top_EC$EC_hi <- rowSums(top_EC[,8:57]> EC_hi)
colnames(top_EC)[59] <- "EC"

# Topsoil organic Carbon
C_lo <- quantile(ref$C/10000, probs=0.25) ## low reference level
C_hi <- quantile(ref$C/10000, probs=0.75) ## high reference level
top_C <- merge(ssid, C, by="SSN")
top_C <- top_C[which(top_C$depth=="top"), ]
top_C <- top_C[!duplicated(top_C[,2]), ]
top_C$C_lo <- rowSums(top_C[,8:57]< C_lo)
top_C$C_hi <- rowSums(top_C[,8:57]> C_hi)
colnames(top_C)[59] <- "C"

# Topsoil Nitrogen
N_lo <- 0.1 ## low reference level
N_hi <- 0.5 ## high reference level
top_N <- merge(ssid, N, by="SSN")
top_N <- top_N[which(top_N$depth=="top"), ]
top_N <- top_N[!duplicated(top_N[,2]), ]
top_N$N_lo <- rowSums(top_N[,8:57]< N_lo)
top_N$N_hi <- rowSums(top_N[,8:57]> N_hi)
colnames(top_N)[59] <- "N"

# Topsoil Boron
B_lo <- 0.5 ## low reference level
B_hi <- 4.0 ## high reference level
top_B <- merge(ssid, B, by="SSN")
top_B <- top_B[which(top_B$depth=="top"), ]
top_B <- top_B[!duplicated(top_B[,2]), ]
top_B$B_lo <- rowSums(top_B[,8:57]< B_lo)
top_B$B_hi <- rowSums(top_B[,8:57]> B_hi)
colnames(top_B)[59] <- "B"

# Topsoil Magnesium
Mg_lo <- quantile(ref$Mg, probs=0.25) ## low reference level
Mg_hi <- quantile(ref$Mg, probs=0.75) ## high reference level
top_Mg <- merge(ssid, Mg, by="SSN")
top_Mg <- top_Mg[which(top_Mg$depth=="top"), ]
top_Mg <- top_Mg[!duplicated(top_Mg[,2]), ]
top_Mg$Mg_lo <- rowSums(top_Mg[,8:57]< Mg_lo)
top_Mg$Mg_hi <- rowSums(top_Mg[,8:57]> Mg_hi)
colnames(top_Mg)[59] <- "Mg"

# Topsoil Phosporus
P_lo <- 15 ## low reference level
P_hi <- 150 ## high reference level
top_P <- merge(ssid, P, by="SSN")
top_P <- top_P[which(top_P$depth=="top"), ]
top_P <- top_P[!duplicated(top_P[,2]), ]
top_P$P_lo <- rowSums(top_P[,8:57]< P_lo)
top_P$P_hi <- rowSums(top_P[,8:57]> P_hi)
colnames(top_P)[59] <- "P"

# Topsoil Sulfur
S_lo <- 10 ## low reference level
S_hi <- 100 ## high reference level
top_S <- merge(ssid, S, by="SSN")
top_S <- top_S[which(top_S$depth=="top"), ]
top_S <- top_S[!duplicated(top_S[,2]), ]
top_S$S_lo <- rowSums(top_S[,8:57]< S_lo)
top_S$S_hi <- rowSums(top_S[,8:57]> S_hi)
colnames(top_S)[59] <- "S"

# Topsoil Potassium
K_lo <- 90 ## low reference level
K_hi <- 900 ## high reference level
top_K <- merge(ssid, K, by="SSN")
top_K <- top_K[which(top_K$depth=="top"), ]
top_K <- top_K[!duplicated(top_K[,2]), ]
top_K$K_lo <- rowSums(top_K[,8:57]< K_lo)
top_K$K_hi <- rowSums(top_K[,8:57]> K_hi)
colnames(top_K)[59] <- "K"

# Topsoil Calcium
Ca_lo <- quantile(ref$Ca, probs=0.25) ## low reference level
Ca_hi <- quantile(ref$Ca, probs=0.75) ## high reference level
top_Ca <- merge(ssid, Ca, by="SSN")
top_Ca <- top_Ca[which(top_Ca$depth=="top"), ]
top_Ca <- top_Ca[!duplicated(top_Ca[,2]), ]
top_Ca$Ca_lo <- rowSums(top_Ca[,8:57]< Ca_lo)
top_Ca$Ca_hi <- rowSums(top_Ca[,8:57]> Ca_hi)
colnames(top_Ca)[59] <- "Ca"

# Topsoil Manganese
Mn_lo <- quantile(ref$Mn, probs=0.25) ## low reference level
Mn_hi <- quantile(ref$Mn, probs=0.75) ## high reference level
top_Mn <- merge(ssid, Mn, by="SSN")
top_Mn <- top_Mn[which(top_Mn$depth=="top"), ]
top_Mn <- top_Mn[!duplicated(top_Mn[,2]), ]
top_Mn$Mn_lo <- rowSums(top_Mn[,8:57]< Mn_lo)
top_Mn$Mn_hi <- rowSums(top_Mn[,8:57]> Mn_hi)
colnames(top_Mn)[59] <- "Mn"

# Topsoil Iron
Fe_lo <- quantile(ref$Fe, probs=0.25) ## low reference level
Fe_hi <- quantile(ref$Fe, probs=0.75) ## high reference level
top_Fe <- merge(ssid, Fe, by="SSN")
top_Fe <- top_Fe[which(top_Fe$depth=="top"), ]
top_Fe <- top_Fe[!duplicated(top_Fe[,2]), ]
top_Fe$Fe_lo <- rowSums(top_Fe[,8:57]< Fe_lo)
top_Fe$Fe_hi <- rowSums(top_Fe[,8:57]> Fe_hi)
colnames(top_Fe)[59] <- "Fe"

# Topsoil Copper
Cu_lo <- quantile(ref$Cu, probs=0.25) ## low reference level
Cu_hi <- quantile(ref$Cu, probs=0.75) ## high reference level
top_Cu <- merge(ssid, Cu, by="SSN")
top_Cu <- top_Cu[which(top_Cu$depth=="top"), ]
top_Cu <- top_Cu[!duplicated(top_Cu[,2]), ]
top_Cu$Cu_lo <- rowSums(top_Cu[,8:57]< Cu_lo)
top_Cu$Cu_hi <- rowSums(top_Cu[,8:57]> Cu_hi)
colnames(top_Cu)[59] <- "Cu"

# Topsoil Zinc
Zn_lo <- 1 ## low reference level
Zn_hi <- 20 ## high reference level
top_Zn <- merge(ssid, Zn, by="SSN")
top_Zn <- top_Zn[which(top_Zn$depth=="top"), ]
top_Zn <- top_Zn[!duplicated(top_Zn[,2]), ]
top_Zn$Zn_lo <- rowSums(top_Zn[,8:57]< Zn_lo)
top_Zn$Zn_hi <- rowSums(top_Zn[,8:57]> Zn_hi)
colnames(top_Zn)[59] <- "Zn"

# Topsoil Aluminum
Al_lo <- quantile(ref$Al, probs=0.25) ## low reference level
Al_hi <- quantile(ref$Al, probs=0.75) ## high reference level
top_Al <- merge(ssid, Al, by="SSN")
top_Al <- top_Al[which(top_Al$depth=="top"), ]
top_Al <- top_Al[!duplicated(top_Al[,2]), ]
top_Al$Al_lo <- rowSums(top_Al[,8:57]< Al_lo)
top_Al$Al_hi <- rowSums(top_Al[,8:57]> Al_hi)
colnames(top_Al)[59] <- "Al"

# Topsoil Sodium
Na_lo <- quantile(ref$Na, probs=0.25) ## low reference level
Na_hi <- quantile(ref$Na, probs=0.75) ## high reference level
top_Na <- merge(ssid, Na, by="SSN")
top_Na <- top_Na[which(top_Na$depth=="top"), ]
top_Na <- top_Na[!duplicated(top_Na[,2]), ]
top_Na$Na_lo <- rowSums(top_Na[,8:57]< Na_lo)
top_Na$Na_hi <- rowSums(top_Na[,8:57]> Na_hi)
colnames(top_Na)[59] <- "Na"

# Topsoil exchangeable acidity (Hp)
Hp_lo <- quantile(ref$Hp, probs=0.25) ## low reference level
Hp_hi <- quantile(ref$Hp, probs=0.75) ## high reference level
top_Hp <- merge(ssid, Hp, by="SSN")
top_Hp <- top_Hp[which(top_Hp$depth=="top"), ]
top_Hp <- top_Hp[!duplicated(top_Hp[,2]), ]
top_Hp$Hp_lo <- rowSums(top_Hp[,8:57]< Hp_lo)
top_Hp$Hp_hi <- rowSums(top_Hp[,8:57]> Hp_hi)
colnames(top_Hp)[59] <- "Hp"

# Topsoil cum distribution plots ------------------------------------------
par(mfrow=c(2,3), mar=c(5,4.5,1,1))
plo <- 0.25
phi <- 0.75

# Topsoil Nitrogen predictions
plot(ecdf(top_N$N), main="", xlab="N (%)", cex.lab=1.4, ylab="CDF", xlim=c(0,0.3), verticals=T, lty=1, lwd=2, col="red", do.points=F, col.01line = NULL)
abline(v=quantile(ref$N/10000, probs=plo), h=plo, lty=1, col="grey")
abline(v=quantile(ref$N/10000, probs=phi), h=phi, lty=1, col="grey")
plot(ecdf(top_N$p05), add=T, verticals=T, lty=2, lwd=2, col="dark grey", do.points=F, col.01line = NULL)
plot(ecdf(top_N$p95), add=T, verticals=T, lty=2, lwd=2, col="dark grey", do.points=F, col.01line = NULL)

# Topsoil Phosphorus predictions
plot(ecdf(top_P$P), main="", xlab="P (ppm)", cex.lab=1.4, ylab="CDF", xlim=c(0,40), verticals=T, lty=1, lwd=2, col="red", do.points=F, col.01line = NULL)
abline(v=quantile(ref$P, probs=plo), h=plo, lty=1, col="grey")
abline(v=quantile(ref$P, probs=phi), h=phi, lty=1, col="grey")
plot(ecdf(top_P$p05), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)
plot(ecdf(top_P$p95), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)

# Topsoil Potassium predictions
plot(ecdf(top_K$K), main="", xlab="K (ppm)", cex.lab=1.4, ylab="CDF", xlim=c(0,300), verticals=T, lty=1, lwd=2, col="red", do.points=F, col.01line = NULL)
abline(v=quantile(ref$K, probs=plo), h=plo, lty=1, col="grey")
abline(v=quantile(ref$K, probs=phi), h=phi, lty=1, col="grey")
plot(ecdf(top_K$p05), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)
plot(ecdf(top_K$p95), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)

# Topsoil Sulfur predictions
plot(ecdf(top_S$S), main="", xlab="S (ppm)", cex.lab=1.4, ylab="CDF", xlim=c(0,40), verticals=T, lty=1, lwd=2, col="red", do.points=F, col.01line = NULL)
abline(v=quantile(ref$S, probs=plo), h=plo, lty=1, col="grey")
abline(v=quantile(ref$S, probs=phi), h=phi, lty=1, col="grey")
plot(ecdf(top_S$p05), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)
plot(ecdf(top_S$p95), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)

# Topsoil Boron predictions
plot(ecdf(top_B$B), main="", xlab="B (ppm)", cex.lab=1.4, ylab="CDF", xlim=c(0,1), verticals=T, lty=1, lwd=2, col="red", do.points=F, col.01line = NULL)
abline(v=quantile(ref$B, probs=plo), h=plo, lty=1, col="grey")
abline(v=quantile(ref$B, probs=phi), h=phi, lty=1, col="grey")
plot(ecdf(top_B$p05), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)
plot(ecdf(top_B$p95), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)

# Topsoil Zinc predictions
plot(ecdf(top_Zn$Zn), main="", xlab="Zn (ppm)", cex.lab=1.4, ylab="CDF", xlim=c(0,4), verticals=T, lty=1, lwd=2, col="red", do.points=F, col.01line = NULL)
abline(v=quantile(ref$Zn, probs=plo), h=plo, lty=1, col="grey")
abline(v=quantile(ref$Zn, probs=phi), h=phi, lty=1, col="grey")
plot(ecdf(top_Zn$p05), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)
plot(ecdf(top_Zn$p95), add=T, verticals=T, lty=1, lwd=2, col="dark grey", do.points=F, col.01line = NULL)

par(mfrow=c(1,1))

# Write file --------------------------------------------------------------
mirpred <- cbind(top_pH[c(1:7,59,61:62)], top_EC[c(59,61:62)], top_C[c(59,61:62)], top_N[c(59,61:62)], top_B[c(59,61:62)],
                 top_Mg[c(59,61:62)], top_P[c(59,61:62)], top_S[c(59,61:62)], top_K[c(59,61:62)], top_Ca[c(59,61:62)], 
                 top_Mn[c(59,61:62)], top_Fe[c(59,61:62)], top_Cu[c(59,61:62)], top_Zn[c(59,61:62)], top_Al[c(59,61:62)],
                 top_Na[c(59,61:62)], top_Hp[c(59,61:62)])

# Project coords to Lambert Azimuthal Equal Area (laea) CRS
mirpred.proj <- as.data.frame(project(cbind(mirpred$lon, mirpred$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(mirpred.proj) <- c("x","y") ## laea coordinates
mirpred <- cbind(mirpred, mirpred.proj)
write.csv(mirpred, "Top_MIR_pred.csv", row.names = F)



