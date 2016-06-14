#' Mehlich-3 (M3) data checks by comparison to the current AfSIS reference data
#' M. Walsh, May 2016

# install.packages(c("downloader","compositions"), dependencies=T)
require(downloader)
require(compositions)

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("M3_data", showWarnings=F)
setwd("./M3_data")

# Download AfSIS reference data
download("https://www.dropbox.com/s/4cfu3crtbcslleq/AfSIS_M3_top.csv.zip?dl=0", "AfSIS_M3_top.csv.zip", mode="wb")
unzip("AfSIS_M3_top.csv.zip", overwrite=T)
ref <- read.table("AfSIS_M3_top.csv", header=T, sep=",") ## AfSIS-M3 topsoil reference data in mg/kg

# Copy M3 sample data into ./M3_data ... working directory
unzip("OCP_M3.zip", overwrite=T)
samp <- read.table("OCP_M3_rerun.csv", header=T, sep=",") ## topsoil sample M3 data in mg/kg

# Compositional data setup
dat <- rbind(ref, samp)
vars <- c("B","Mg","P","S","K","Ca","Mn","Fe","Cu","Zn")
qdat <- na.omit(dat[vars])
qdat$Fv <- 1000000-rowSums(qdat[vars]) ## calculates "fill value" (Fv), in mg/kg soil

# Centered log ratio (clr) transform
cvar <- c("B","Mg","P","S","K","Ca","Mn","Fe","Cu","Zn","Fv")
cdat <- acomp(qdat[cvar])
clrt <- as.data.frame(clr(cdat)) ## clr transform
cdat <- cbind(dat[1:2], clrt) ## generates compositional data frame seperated by reference (RS) values

# GLM tests ---------------------------------------------------------------
ftest <- glm(RS ~ Fv, family=binomial(link="logit"), cdat) ## differences in fill values (Fv)
summary(ftest)
cdat$fOR <- predict(ftest, type="link", cdat) ## log odds ratio of the difference between reference and sample measurements

ctest <- glm(RS ~ B+Mg+P+S+K+Ca+Mn+Fe+Cu+Zn, family=binomial(link="logit"), cdat)
summary(ctest)
cdat$cOR <- predict(ctest, type="link", cdat) ## log odds ratio of the difference between reference and sample measurements

# Diagnostic plots
par(mfrow=c(1,2))
boxplot(fOR~RS, notch=T, ylab="Fill value index", cdat)
boxplot(cOR~RS, notch=T, ylab="Compositional index", cdat)
par(mfrow=c(1,1))
plot(cOR~fOR, col=ifelse(RS=="R", 4, 2), cex=0.6, xlab="Fill value index", ylab="Compositional index", cdat)
