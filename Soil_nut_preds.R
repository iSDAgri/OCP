#' Soil nutrient map predictions from soil MIR and gridded data
#' M. Walsh, December 2016

# Required packages
# install.packages(c("devtools","caret","plyr","doParallel" ...)), dependencies=TRUE)
suppressPackageStartupMessages({
  require(devtools)
  require(caret)
  require(plyr)
  require(doParallel)
  require(raster)
  require(randomForest)
  require(gbm)
  require(nnet)
  require(bartMachine)
  require(glmnet)
})

# Data setup --------------------------------------------------------------
# SourceURL <- "https://raw.githubusercontent.com/mgwalsh/OCP/blob/master/Spatial_pred_setup.R"
# source_url(SourceURL)

# Nutrient reference levels (in ppm)
rN <- 1500
rP <- 30
rK <- 190
rS <- 20
rB <- 0.8
rZn <- 1.5

# Oxide equivalents
oP <- 2.2913 ## P to P2O5
oK <- 1.2046 ## K to K2O

# Nutrient mass estimates (kg/ha) in delta notation relative to specified reference levels
nut$dN <- 2*nut$BD20*nut$N*10000-2*nut$BD20*rN ## convert from % to ppm
nut$dP <- 2*nut$BD20*nut$P*oP-2*nut$BD20*rP ## in P2O5 oxide equivalents
nut$dK <- 2*nut$BD20*nut$K*oK-2*nut$BD20*rK ## in K2O oxide equivalents
nut$dS <- 2*nut$BD20*nut$S-2*nut$BD20*rS
nut$dB <- 2*nut$BD20*nut$B-2*nut$BD20*rB
nut$dZn <- 2*nut$BD20*nut$Zn-2*nut$BD20*rZn

# Variable setup
y <- nut$dP ## set nutrient to be predicted here
x <- nut[,11:50] ## set features here

# bartMachine model -------------------------------------------------------
# options(java.parameters = "-Xmx2500M")
# set_bart_machine_num_cores(4)

# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Model setup
set.seed(1385321)
tc <- trainControl(method = "cv", allowParallel = T)
bar <- train(x, y,
             method = "bartMachine", 
             preProc = c("center", "scale"), 
             trControl = tc,
             tuneLength = 2,
             verbose = FALSE,
             seed = 1)
print(bar)
bar_pred <- predict(grids, bar)
rm("bar")

stopCluster(mc)

# RF model ----------------------------------------------------------------
# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Model setup
set.seed(1385321)
tc <- trainControl(method = "cv")
tg <- expand.grid(mtry=seq(2, 10, by=1))
rfo <- train(x, y,
             preProc = c("center", "scale"),
             method = "rf",
             ntree = 501,
             tuneGrid = tg,
             trControl = tc)
print(rfo)
rfo_pred <- predict(grids, rfo) ## predict map
rm("rfo")

stopCluster(mc)

# GBM model ---------------------------------------------------------------
# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Model setup
set.seed(1385321)
tc <- trainControl(method = "repeatedcv", repeats=5)
tg <- expand.grid(.n.trees=seq(10, 200, by=5), 
                  .interaction.depth = 5,
                  .shrinkage = 0.1,
                  .n.minobsinnode = 10)
gbm <- train(x, y, 
             method = "gbm", 
             preProc = c("center", "scale"),
             trControl = tc,
             tuneGrid = tg)
print(gbm)
gbm_pred <- predict(grids, gbm) ## predict map
rm("gbm")

stopCluster(mc)

# NNET models -------------------------------------------------------------
# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Model setup
set.seed(1385321)
tc <- trainControl(method = "cv")
tg <- expand.grid(layer1 = 2:10,
                  layer2 = 0,
                  layer3 = 0,
                  hidden_dropout = 0,
                  visible_dropout = 0)
net <- train(x, y, 
             method = "dnn", 
             preProc = c("center", "scale"), 
             trControl = tc,
             tuneGrid = tg)
print(net)
net_pred <- predict(grids, net) ## predict map
rm("net")

stopCluster(mc)

# Model stacking setup ----------------------------------------------------
pwetv <- as.data.frame(cbind(rfo_wet, gbm_wet, dnn_wet, bar_wet))
names(pwetv) <- c("L", "RFO", "GBM", "NNET", "BART")

# Write data files --------------------------------------------------------
write.csv(pwetv, "pwetv.csv", row.names=F)
write.csv(pmirv, "pmirv.csv", row.names=F)

# Remove extraneous objects from memory -----------------------------------
# rm(list=setdiff(ls(), c("pwetv", "pmirv")))

# Model stacking ----------------------------------------------------------
# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Model setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = TRUE, summaryFunction = twoClassSummary,
                   allowParallel = T)
wet.ens <- train(HL ~ ., data = pwetv,
                 method = "glmnet",
                 family = "binomial",
                 metric = "ROC",
                 trControl = tc)
print(wet.ens)
wet.imp <- varImp(wet.ens)
plot(wet.imp, top=4, col="black", cex=1.2, xlab="Model importance in ensemble prediction")
ens_wet <- predict(wet.ens, pwetv, type = "prob")

stopCluster(mc)
