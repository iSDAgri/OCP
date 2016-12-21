#' Central Nigeria cropland and maize sentinel site survey summaries 
#' M. Walsh, September 2016

# Required packages
# install.packages(c("downloader","rgdal","proj4","leaderCluster","arm")), dependencies=TRUE)
require(downloader)
require(rgdal)
require(proj4)
require(leaderCluster)
require(arm)

# Data downloads -----------------------------------------------------------
# Create a "Data" folder in your current working directory
dir.create("MS_data", showWarnings=F)
setwd("./MS_data")

# download GeoSurvey data
download("https://www.dropbox.com/s/6gr2s8s7zgj7hr7/GS_CNG_2016.csv.zip?dl=0", "GS_CNG_2016.csv.zip", mode="wb")
unzip("GS_CNG_2016.csv.zip", overwrite=T)
gsdat <- read.table("GS_CNG_2016.csv", header=T, sep=",")

# download MobileSurvey data
download("https://www.dropbox.com/s/7gxvhg0ztr59nqw/OCP_maize.csv.zip?dl=0", "OCP_maize.csv.zip", mode="wb")
unzip("OCP_maize.csv.zip", overwrite=T)
msdat <- read.table("OCP_maize.csv", header=T, sep=",")

# GeoSurvey area estimates ------------------------------------------------
roi <- 239504 # total ROI area (in km^2)
gsdat <- gsdat[which(gsdat$Sample==1),]
n <- nrow(gsdat) # total number of GeoSurvey observations

# cropland presence
p <- nrow(gsdat[which(gsdat$CRP=="Y"),]) # number of observations with cropland presence
cp <- p/n # proportion
se <- sqrt(cp*(1-cp)/n) # standard error
cp.95 <- cp + qnorm(c(0.025, 0.5, 0.975))*se # 95% confidence interval
cp.roi <- roi * cp.95 * 1/10000

# presence of buildings
p <- nrow(gsdat[which(gsdat$RSP=="Y"),]) # number of observations with presence of buildings
bp <- p/n # proportion
se <- sqrt(bp*(1-bp)/n) # standard error
bp.95 <- bp + qnorm(c(0.025, 0.5, 0.975)) * se # 95% confidence interval
bp.roi <- roi * bp.95 * 1/10000

# presence of of woody vegetation cover >60%
p <- nrow(gsdat[which(gsdat$WCP=="Y"),]) # number of observations with woody cover presence
wp <- p/n # proportion
se <- sqrt(wp*(1-wp)/n) # standard error
wp.95 <- wp + qnorm(c(0.025, 0.5, 0.975)) * se # 95% confidence interval
wp.est <- roi * wp.95 * 1/10000

# Generate coordinate reference and Site IDs ------------------------------
# Project MobileSurvey coords to Africa LAEA from LonLat
msdat.laea <- as.data.frame(project(cbind(msdat$lon, msdat$lat), "+proj=laea +ellps=WGS84 +lon_0=20 +lat_0=5 +units=m +no_defs"))
colnames(msdat.laea) <- c("x","y")
msdat <- cbind(msdat, msdat.laea)
loc <- as.matrix(msdat[,8:9])
sid <- leaderCluster(loc, radius=12000, distance="L2")$cluster_id
msdat <- cbind(msdat, sid)
ssdat <- aggregate(.~sid, data = msdat, mean)

# Plot of raw sentinel site data ------------------------------------------
par(mfrow=c(1,2), mar=c(5,4.5,1,1))
plot(crop~gs_cp, ssdat, xlab="Proportion GS cropland", ylab="Proportion MS cropland", xlim=c(0,1), ylim=c(0,1))
abline(c(0,1), lwd=1)
plot(maize~crop, ssdat, xlab="Proportion MS cropland", ylab="Proportion MS Maize", xlim=c(0,1), ylim=c(0,1))
abline(c(0,1), lwd=1)
par(mfrow=c(1,1))

# GLM models --------------------------------------------------------------
# unconditional presence/absence (p/a) of cropland
c1.glm <- glm(crop~1, family=binomial(link="logit"), data=msdat)
display(c1.glm)

# p/a of cropland conditioned on GeoSurvey observations
c2.glm <- glm(crop~gs_cp+gs_bp+gs_wp, family=binomial(link="logit"), data=msdat)
display(c2.glm)

# unconditional p/a of Maize
m1.glm <- glm(maize~1, family=binomial(link="logit"), data=msdat)
display(m1.glm)

# p/a of Maize conditioned on GeoSurvey observations
m2.glm <- glm(maize~gs_cp+gs_bp+gs_wp, family=binomial(link="logit"), data=msdat)
display(m2.glm)

# GLMER models ------------------------------------------------------------
# p/a of cropland by SID
c1.glmer <- glmer(crop~1+(1|sid), family=binomial(link="logit"), data=msdat)
display(c1.glmer)

# p/a of cropland conditioned on GeoSurvey observations & SID
c2.glmer <- glmer(crop~gs_cp+gs_bp+gs_wp+(1|sid), family=binomial(link="logit"), data=msdat)
display(c2.glmer)
c2.sim <- sim(c2.glmer, n.sims = 1000)
c2.fix <- as.data.frame(fixef(c2.sim))
colnames(c2.fix) <- c("int","cp","bp","wp")
c2.est <- c2.fix$cp*cp+c2.fix$bp*bp+c2.fix$wp*wp+c2.fix$int
c2.q <- quantile(c2.est, probs = c(0.025,0.5,0.975))
c2.est <- exp(c2.q)/(1+exp(c2.q))
c2.est
c2.roi <- (roi * c2.est) * 1/10000 # estimate in Mha

# p/a of Maize by SID
m1.glmer <- glmer(maize~1+(1|sid), family=binomial(link="logit"), data=msdat)
display(m1.glmer)

# p/a of Maize conditioned on GeoSurvey observations & SID
m2.glmer <- glmer(maize~gs_cp+gs_bp+gs_wp+(1|sid), family=binomial(link="logit"), data=msdat)
display(m2.glmer)
m2.sim <- sim(m2.glmer, n.sims = 1000)
m2.fix <- as.data.frame(fixef(m2.sim))
colnames(m2.fix) <- c("int","cp","bp","wp")
m2.est <- m2.fix$cp*cp+m2.fix$bp*bp+m2.fix$wp*wp+m2.fix$int
m2.q <- quantile(m2.est, probs = c(0.025,0.5,0.975))
m2.est <- exp(m2.q)/(1+exp(m2.q))
m2.est
m2.roi <- roi * m2.est * 1/10000 # estimate in Mha

# p/a of Sorghum conditioned on GeoSurvey observations & SID
s2.glmer <- glmer(sorgh~gs_cp+gs_bp+gs_wp+(1|sid), family=binomial(link="logit"), data=msdat)
display(s2.glmer)
s2.sim <- sim(s2.glmer, n.sims = 1000)
s2.fix <- as.data.frame(fixef(s2.sim))
colnames(s2.fix) <- c("int","cp","bp","wp")
s2.est <- s2.fix$cp*cp+s2.fix$bp*bp+s2.fix$wp*wp+s2.fix$int
s2.q <- quantile(s2.est, probs = c(0.025,0.5,0.975))
s2.est <- exp(s2.q)/(1+exp(s2.q))
s2.est
s2.roi <- roi * s2.est * 1/10000 # estimate in Mha
