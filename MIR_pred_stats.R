#' Data summaries of MIR-M3 prediction results
#' M. Walsh, June 2016

# install.packages(c("downloader","compositions"), dependencies=T)
require(downloader)

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("MIR_data", showWarnings=F)
setwd("./MIR_data")

# Download MIR-M3 prediction data
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

# Sample ID codebook
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
# Topsoil pH predictions
top_pH <- merge(ssid, pH, by="SSN")
top_pH <- top_pH[which(top_pH$depth=="top"), ]
top_pH <- top_pH[!duplicated(top_pH[,2]), ]

# Topsoil electrical conductivity (EC) predictions
top_EC <- merge(ssid, EC, by="SSN")
top_EC <- top_EC[which(top_EC$depth=="top"), ]
top_EC <- top_EC[!duplicated(top_EC[,2]), ]

top_C <- merge(ssid, C, by="SSN")
top_C <- top_C[which(top_C$depth=="top"), ]
top_C <- top_C[!duplicated(top_C[,2]), ]

# Topsoil Nitrogen predictions
top_N <- merge(ssid, N, by="SSN")
top_N <- top_N[which(top_N$depth=="top"), ]
top_N <- top_N[!duplicated(top_N[,2]), ]

# Topsoil Boron predictions
top_B <- merge(ssid, B, by="SSN")
top_B <- top_B[which(top_B$depth=="top"), ]
top_B <- top_B[!duplicated(top_B[,2]), ]

# Topsoil Magnesium predictions
top_Mg <- merge(ssid, Mg, by="SSN")
top_Mg <- top_Mg[which(top_Mg$depth=="top"), ]
top_Mg <- top_Mg[!duplicated(top_Mg[,2]), ]

# Topsoil Phosporus predictions
top_P <- merge(ssid, P, by="SSN")
top_P <- top_P[which(top_P$depth=="top"), ]
top_P <- top_P[!duplicated(top_P[,2]), ]

# Topsoil Sulfur predictions
top_S <- merge(ssid, S, by="SSN")
top_S <- top_S[which(top_S$depth=="top"), ]
top_S <- top_S[!duplicated(top_S[,2]), ]

# Topsoil Potassium predictions
top_K <- merge(ssid, K, by="SSN")
top_K <- top_K[which(top_K$depth=="top"), ]
top_K <- top_K[!duplicated(top_K[,2]), ]

# Topsoil Calcium predictions
top_Ca <- merge(ssid, Ca, by="SSN")
top_Ca <- top_Ca[which(top_Ca$depth=="top"), ]
top_Ca <- top_Ca[!duplicated(top_Ca[,2]), ]

# Topsoil Manganese predictions
top_Mn <- merge(ssid, Mn, by="SSN")
top_Mn <- top_Mn[which(top_Mn$depth=="top"), ]
top_Mn <- top_Mn[!duplicated(top_Mn[,2]), ]

# Topsoil Iron predictions
top_Fe <- merge(ssid, Fe, by="SSN")
top_Fe <- top_Fe[which(top_Fe$depth=="top"), ]
top_Fe <- top_Fe[!duplicated(top_Fe[,2]), ]

# Topsoil Copper predictions
top_Cu <- merge(ssid, Cu, by="SSN")
top_Cu <- top_Cu[which(top_Cu$depth=="top"), ]
top_Cu <- top_Cu[!duplicated(top_Cu[,2]), ]

# Topsoil Zinc predictions
top_Zn <- merge(ssid, Zn, by="SSN")
top_Zn <- top_Zn[which(top_Zn$depth=="top"), ]
top_Zn <- top_Zn[!duplicated(top_Zn[,2]), ]

# Topsoil Aluminum predictions
top_Al <- merge(ssid, Al, by="SSN")
top_Al <- top_Al[which(top_Al$depth=="top"), ]
top_Al <- top_Al[!duplicated(top_Al[,2]), ]

# Topsoil Sodium predictions
top_Na <- merge(ssid, Na, by="SSN")
top_Na <- top_Na[which(top_Na$depth=="top"), ]
top_Na <- top_Na[!duplicated(top_Na[,2]), ]

# Topsoil exchangeable acidity (Hp) predictions
top_Hp <- merge(ssid, Hp, by="SSN")
top_Hp <- top_Hp[which(top_Hp$depth=="top"), ]
top_Hp <- top_Hp[!duplicated(top_Hp[,2]), ]

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
# Topsoil pH predictions
names(top_pH)[58:60] <- c("p.5", "p.50", "p.95")
top_pH <- reshape(top_pH, direction="long", varying=58:60, idvar="ssid", v.names="pH", timevar="plevel") ## long format
top_pH$lo <- ifelse(top_pH$pH < quantile(ref$pH, probs=plo), 1, 0) ## identifies low levels
top_pH$hi <- ifelse(top_pH$pH > quantile(ref$pH, probs=phi), 1, 0) ## identifies high levels
top_pH <- top_pH[,c(1:7,58:61)]

# Topsoil electrical conductivity (EC) predictions
names(top_EC)[58:60] <- c("p.5", "p.50", "p.95")
top_EC <- reshape(top_EC, direction="long", varying=58:60, idvar="ssid", v.names="EC", timevar="plevel") ## long format
top_EC$lo <- ifelse(top_EC$EC < quantile(ref$EC*1000, probs=plo), 1, 0) ## identifies low levels
top_EC$hi <- ifelse(top_EC$EC > quantile(ref$EC*1000, probs=phi), 1, 0) ## identifies high levels
top_EC <- top_EC[,c(1:7,58:61)]

# Topsoil organic Carbon predictions
names(top_C)[58:60] <- c("p.5", "p.50", "p.95")
top_C <- reshape(top_C, direction="long", varying=58:60, idvar="ssid", v.names="C", timevar="plevel") ## long format
top_C$lo <- ifelse(top_C$C < quantile(ref$C/10000, probs=plo), 1, 0) ## identifies low levels
top_C$hi <- ifelse(top_C$C > quantile(ref$C/10000, probs=phi), 1, 0) ## identifies high levels
top_C <- top_C[,c(1:7,58:61)]

names(top_N)[58:60] <- c("p.5", "p.50", "p.95")
top_N <- reshape(top_N, direction="long", varying=58:60, idvar="ssid", v.names="N", timevar="plevel") ## long format
top_N$lo <- ifelse(top_N$N < quantile(ref$N/10000, probs=plo), 1, 0) ## identifies low levels
top_N$hi <- ifelse(top_N$N > quantile(ref$N/10000, probs=phi), 1, 0) ## identifies high levels
top_N <- top_N[,c(1:7,58:61)]

# Topsoil Boron predictions
names(top_B)[58:60] <- c("p.5", "p.50", "p.95")
top_B <- reshape(top_B, direction="long", varying=58:60, idvar="ssid", v.names="B", timevar="plevel") ## long format
top_B$lo <- ifelse(top_B$B < quantile(ref$B, probs=plo), 1, 0) ## identifies low levels
top_B$hi <- ifelse(top_B$B > quantile(ref$B, probs=phi), 1, 0) ## identifies high levels
top_B <- top_B[,c(1:7,58:61)]

# Topsoil Magnesium predictions
names(top_Mg)[58:60] <- c("p.5", "p.50", "p.95")
top_Mg <- reshape(top_Mg, direction="long", varying=58:60, idvar="ssid", v.names="Mg", timevar="plevel") ## long format
top_Mg$lo <- ifelse(top_Mg$Mg < quantile(ref$Mg, probs=plo), 1, 0) ## identifies low levels
top_Mg$hi <- ifelse(top_Mg$Mg > quantile(ref$Mg, probs=phi), 1, 0) ## identifies high levels
top_Mg <- top_Mg[,c(1:7,58:61)]

# Topsoil Phosphorus predictions
names(top_P)[58:60] <- c("p.5", "p.50", "p.95")
top_P <- reshape(top_P, direction="long", varying=58:60, idvar="ssid", v.names="P", timevar="plevel") ## long format
top_P$lo <- ifelse(top_P$P < quantile(ref$P, probs=plo), 1, 0) ## identifies low levels
top_P$hi <- ifelse(top_P$P > quantile(ref$P, probs=phi), 1, 0) ## identifies high levels
top_P <- top_P[,c(1:7,58:61)]

# Topsoil Sulfur predictions
names(top_S)[58:60] <- c("p.5", "p.50", "p.95")
top_S <- reshape(top_S, direction="long", varying=58:60, idvar="ssid", v.names="S", timevar="plevel") ## long format
top_S$lo <- ifelse(top_S$S < quantile(ref$S, probs=plo), 1, 0) ## identifies low levels
top_S$hi <- ifelse(top_S$S > quantile(ref$S, probs=phi), 1, 0) ## identifies high levels
top_S <- top_S[,c(1:7,58:61)]

# Topsoil Potassium predictions
names(top_K)[58:60] <- c("p.5", "p.50", "p.95")
top_K <- reshape(top_K, direction="long", varying=58:60, idvar="ssid", v.names="K", timevar="plevel") ## long format
top_K$lo <- ifelse(top_K$K < quantile(ref$K, probs=plo), 1, 0) ## identifies low levels
top_K$hi <- ifelse(top_K$K > quantile(ref$K, probs=phi), 1, 0) ## identifies high levels
top_K <- top_K[,c(1:7,58:61)]

# Topsoil Calcium predictions
names(top_Ca)[58:60] <- c("p.5", "p.50", "p.95")
top_Ca <- reshape(top_Ca, direction="long", varying=58:60, idvar="ssid", v.names="Ca", timevar="plevel") ## long format
top_Ca$lo <- ifelse(top_Ca$Ca < quantile(ref$Ca, probs=plo), 1, 0) ## identifies low levels
top_Ca$hi <- ifelse(top_Ca$Ca > quantile(ref$Ca, probs=phi), 1, 0) ## identifies high levels
top_Ca <- top_Ca[,c(1:7,58:61)]

# Topsoil Manganese predictions
names(top_Mn)[58:60] <- c("p.5", "p.50", "p.95")
top_Mn <- reshape(top_Mn, direction="long", varying=58:60, idvar="ssid", v.names="Mn", timevar="plevel") ## long format
top_Mn$lo <- ifelse(top_Mn$Mn < quantile(ref$Mn, probs=plo), 1, 0) ## identifies low levels
top_Mn$hi <- ifelse(top_Mn$Mn > quantile(ref$Mn, probs=phi), 1, 0) ## identifies high levels
top_Mn <- top_Mn[,c(1:7,58:61)]

# Topsoil Iron predictions
names(top_Fe)[58:60] <- c("p.5", "p.50", "p.95")
top_Fe <- reshape(top_Fe, direction="long", varying=58:60, idvar="ssid", v.names="Fe", timevar="plevel") ## long format
top_Fe$lo <- ifelse(top_Fe$Fe < quantile(ref$Fe, probs=plo), 1, 0) ## identifies low levels
top_Fe$hi <- ifelse(top_Fe$Fe > quantile(ref$Fe, probs=phi), 1, 0) ## identifies high levels
top_Fe <- top_Fe[,c(1:7,58:61)]

# Topsoil Copper predictions
names(top_Cu)[58:60] <- c("p.5", "p.50", "p.95")
top_Cu <- reshape(top_Cu, direction="long", varying=58:60, idvar="ssid", v.names="Cu", timevar="plevel") ## long format
top_Cu$lo <- ifelse(top_Cu$Cu < quantile(ref$Cu, probs=plo), 1, 0) ## identifies low levels
top_Cu$hi <- ifelse(top_Cu$Cu > quantile(ref$Cu, probs=phi), 1, 0) ## identifies high levels
top_Cu <- top_Cu[,c(1:7,58:61)]

# Topsoil Zinc predictions
names(top_Zn)[58:60] <- c("p.5", "p.50", "p.95")
top_Zn <- reshape(top_Zn, direction="long", varying=58:60, idvar="ssid", v.names="Zn", timevar="plevel") ## long format
top_Zn$lo <- ifelse(top_Zn$Zn < quantile(ref$Zn, probs=plo), 1, 0) ## identifies low levels
top_Zn$hi <- ifelse(top_Zn$Zn > quantile(ref$Zn, probs=phi), 1, 0) ## identifies high levels
top_Zn <- top_Zn[,c(1:7,58:61)]

# Topsoil Aluminum predictions
names(top_Al)[58:60] <- c("p.5", "p.50", "p.95")
top_Al <- reshape(top_Al, direction="long", varying=58:60, idvar="ssid", v.names="Al", timevar="plevel") ## long format
top_Al$lo <- ifelse(top_Al$Al < quantile(ref$Al, probs=plo), 1, 0) ## identifies low levels
top_Al$hi <- ifelse(top_Al$Al > quantile(ref$Al, probs=phi), 1, 0) ## identifies high levels
top_Al <- top_Al[,c(1:7,58:61)]

# Topsoil Sodium predictions
names(top_Na)[58:60] <- c("p.5", "p.50", "p.95")
top_Na <- reshape(top_Na, direction="long", varying=58:60, idvar="ssid", v.names="Na", timevar="plevel") ## long format
top_Na$lo <- ifelse(top_Na$Na < quantile(ref$Na, probs=plo), 1, 0) ## identifies low levels
top_Na$hi <- ifelse(top_Na$Na > quantile(ref$Na, probs=phi), 1, 0) ## identifies high levels
top_Na <- top_Na[,c(1:7,58:61)]

# Topsoil exchangeable acidity (Hp) predictions
names(top_Hp)[58:60] <- c("p.5", "p.50", "p.95")
top_Hp <- reshape(top_Hp, direction="long", varying=58:60, idvar="ssid", v.names="Hp", timevar="plevel") ## long format
top_Hp$lo <- ifelse(top_Hp$Hp < quantile(ref$Hp, probs=plo), 1, 0) ## identifies low levels
top_Hp$hi <- ifelse(top_Hp$Hp > quantile(ref$Hp, probs=phi), 1, 0) ## identifies high levels
top_Hp <- top_Hp[,c(1:7,58:61)]

# Hi/Lo summaries ---------------------------------------------------------
# Topsoil pH
pH_lo <- table(top_pH$plevel, top_pH$lo)
prop.table(pH_lo,1)*100
pH_hi <- table(top_pH$plevel, top_pH$hi)
prop.table(pH_hi,1)*100

# Topsoil electrical conductivity (EC)
EC_lo <- table(top_EC$plevel, top_EC$lo)
prop.table(EC_lo,1)*100
EC_hi <- table(top_EC$plevel, top_EC$hi)
prop.table(EC_hi,1)*100

# Topsoil organic Carbon
C_lo <- table(top_C$plevel, top_C$lo)
prop.table(C_lo,1)*100
C_hi <- table(top_C$plevel, top_C$hi)
prop.table(C_hi,1)*100

# Topsoil Nitrogen
N_lo <- table(top_N$plevel, top_N$lo)
prop.table(N_lo,1)*100
N_hi <- table(top_N$plevel, top_N$hi)
prop.table(N_hi,1)*100

# Topsoil Boron
B_lo <- table(top_B$plevel, top_B$lo)
prop.table(B_lo,1)*100
B_hi <- table(top_B$plevel, top_B$hi)
prop.table(B_hi,1)*100

# Topsoil Magnesium
Mg_lo <- table(top_Mg$plevel, top_Mg$lo)
prop.table(Mg_lo,1)*100
Mg_hi <- table(top_Mg$plevel, top_Mg$hi)
prop.table(Mg_hi,1)*100

# Topsoil Phosphorus
P_lo <- table(top_P$plevel, top_P$lo)
prop.table(P_lo,1)*100
P_hi <- table(top_P$plevel, top_P$hi)
prop.table(P_hi,1)*100

# Topsoil Sulfur
S_lo <- table(top_S$plevel, top_S$lo)
prop.table(S_lo,1)*100
S_hi <- table(top_S$plevel, top_S$hi)
prop.table(S_hi,1)*100

# Topsoil Potassium
K_lo <- table(top_K$plevel, top_K$lo)
prop.table(K_lo,1)*100
K_hi <- table(top_K$plevel, top_K$hi)
prop.table(K_hi,1)*100

# Topsoil Calcium
Ca_lo <- table(top_Ca$plevel, top_Ca$lo)
prop.table(Ca_lo,1)*100
Ca_hi <- table(top_Ca$plevel, top_Ca$hi)
prop.table(Ca_hi,1)*100

# Topsoil Manganese
Mn_lo <- table(top_Mn$plevel, top_Mn$lo)
prop.table(Mn_lo,1)*100
Mn_hi <- table(top_Mn$plevel, top_Mn$hi)
prop.table(Mn_hi,1)*100

# Topsoil Iron
Fe_lo <- table(top_Fe$plevel, top_Fe$lo)
prop.table(Fe_lo,1)*100
Fe_hi <- table(top_Fe$plevel, top_Fe$hi)
prop.table(Fe_hi,1)*100

# Topsoil Copper
Cu_lo <- table(top_Cu$plevel, top_Cu$lo)
prop.table(Cu_lo,1)*100
Cu_hi <- table(top_Cu$plevel, top_Cu$hi)
prop.table(Cu_hi,1)*100

# Topsoil Zinc
Zn_lo <- table(top_Zn$plevel, top_Zn$lo)
prop.table(Zn_lo,1)*100
Zn_hi <- table(top_Zn$plevel, top_Zn$hi)
prop.table(Zn_hi,1)*100

# Topsoil Aluminum
Al_lo <- table(top_Al$plevel, top_Al$lo)
prop.table(Al_lo,1)*100
Al_hi <- table(top_Al$plevel, top_Al$hi)
prop.table(Al_hi,1)*100

# Topsoil Sodium
Na_lo <- table(top_Na$plevel, top_Na$lo)
prop.table(Na_lo,1)*100
Na_hi <- table(top_Na$plevel, top_Na$hi)
prop.table(Na_hi,1)*100

# Topsoil exchangeable acidity (Hp)
Hp_lo <- table(top_Hp$plevel, top_Hp$lo)
prop.table(Hp_lo,1)*100
Hp_hi <- table(top_Hp$plevel, top_Hp$hi)
prop.table(Hp_hi,1)*100

# Write long MIR prediction file ------------------------------------------
lpred <- cbind(top_pH[,1:9], top_EC["EC"], top_C["C"], top_N["N"], top_B["B"], top_Mg["Mg"], top_P["P"],
               top_S["S"], top_K["K"], top_Ca["Ca"], top_Mn["Mn"], top_Fe["Fe"], top_Cu["Cu"], top_Zn["Zn"],
               top_Al["Al"], top_Na["Na"], top_Hp["Hp"])
write.csv(lpred, "top_MIR_pred_long.csv", row.names=F)

