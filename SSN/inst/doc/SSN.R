## ----results='hide',echo=FALSE, cache=FALSE------------------------------
## Code to run in the background
set.seed(210112)

## ------------------------------------------------------------------------
library("SSN")

## ------------------------------------------------------------------------
file.copy(system.file("lsndata/MiddleFork04.ssn", package = "SSN"),
  to = tempdir(), recursive = TRUE, copy.mode = FALSE)
setwd(tempdir())

## ----eval=FALSE----------------------------------------------------------
## importSSN(Path, predpts = NULL, o.write = FALSE)

## ------------------------------------------------------------------------
mf04p <- importSSN("./MiddleFork04.ssn",
   predpts = "pred1km")

## ------------------------------------------------------------------------
mf04p <- importPredpts(mf04p, "Knapp", "ssn")
mf04p <- importPredpts(mf04p, "CapeHorn", "ssn")

## ----eval=FALSE----------------------------------------------------------
## additive.function(mf04p, VarName, afvName)

## ------------------------------------------------------------------------
names(mf04p@data)

## ------------------------------------------------------------------------
head(mf04p@data[, c("h2oAreaKm2", "afvArea")])

## ------------------------------------------------------------------------
mf04p <- additive.function(mf04p, "h2oAreaKm2", "computed.afv")

## ------------------------------------------------------------------------
names(mf04p@data)

## ------------------------------------------------------------------------
head(mf04p@data[, c("h2oAreaKm2",
   "afvArea", "computed.afv")])
head(getSSNdata.frame(mf04p)[, c("afvArea", "computed.afv")])

## ------------------------------------------------------------------------
createDistMat(mf04p, predpts = "Knapp", o.write = TRUE,
	amongpreds = TRUE)
createDistMat(mf04p, predpts = "CapeHorn", o.write = TRUE,
	amongpreds = TRUE)

## ------------------------------------------------------------------------
distObs <- getStreamDistMat(mf04p)
str(distObs)
distObs$dist.net1[1:5,1:5]

## ------------------------------------------------------------------------
strDistNet2 <- distObs$dist.net2 + t(distObs$dist.net2)
strDistNet2[5:10,5:10]

## ------------------------------------------------------------------------
distPred1km <- getStreamDistMat(mf04p, Name = "pred1km")
str(distPred1km)
distPred1km$dist.net1.a[1:5,1:5]

## ------------------------------------------------------------------------
createDistMat(mf04p, predpts = "CapeHorn", o.write = TRUE, 
  amongpreds = TRUE)
distCape <- getStreamDistMat(mf04p, Name = "CapeHorn")
str(distCape)
distCape$dist.net2[1:5,1:5]

## ------------------------------------------------------------------------
names(mf04p)

## ----LoadData,include=FALSE,results='hide',cache=FALSE-------------------
plot(mf04p, lwdLineCol = "afvArea", lwdLineEx = 10, lineCol = "blue",
   pch = 19, xlab = "x-coordinate (m)", ylab = "y-coordinate (m)",
   asp = 1)

## ----eval=FALSE----------------------------------------------------------
## plot(mf04p, lwdLineCol = "afvArea", lwdLineEx = 10, lineCol = "blue",
##    pch = 19, xlab = "x-coordinate (m)", ylab = "y-coordinate (m)",
##    asp = 1)

## ----plotSpatialStreamNetwork,include=FALSE,results='hide', cache=FALSE----
brks <- plot(mf04p, "Summer_mn", lwdLineCol = "afvArea",
   lwdLineEx = 15, lineCol = "black", xlab =  "x-coordinate" ,
   ylab =  "y-coordinate", asp=1 )

## ----eval=FALSE----------------------------------------------------------
## brks <- plot(mf04p, "Summer_mn", lwdLineCol = "afvArea",
##    lwdLineEx = 15, lineCol = "black", xlab =  "x-coordinate" ,
##    ylab =  "y-coordinate", asp=1 )

## ----plotUsingSP, include=FALSE, results='hide', cache=FALSE-------------
#plot the stream lines
plot(as.SpatialLines(mf04p), col = "blue")
# add the observed locations with size proportional 
# to mean summer temperature
plot(as.SpatialPoints(mf04p), pch = 19, 
  cex = as.SpatialPointsDataFrame(mf04p)$Summer_mn/9 , add = TRUE)
# add the prediction locations on the 1 km spacing
plot(as.SpatialPoints(mf04p, data = "pred1km"), cex = 1.5, add = TRUE)
# add the dense set of points for block prediction on Knapp segment
plot(as.SpatialPoints(mf04p, data = "Knapp"), pch = 19, cex = 0.3, 
  col = "red", add = TRUE)

## ----eval=FALSE----------------------------------------------------------
## #plot the stream lines
## plot(as.SpatialLines(mf04p), col = "blue")
## # add the observed locations with size proportional
## # to mean summer temperature
## plot(as.SpatialPoints(mf04p), pch = 19,
##   cex = as.SpatialPointsDataFrame(mf04p)$Summer_mn/9 , add = TRUE)
## # add the prediction locations on the 1 km spacing
## plot(as.SpatialPoints(mf04p, data = "pred1km"), cex = 1.5, add = TRUE)
## # add the dense set of points for block prediction on Knapp segment
## plot(as.SpatialPoints(mf04p, data = "Knapp"), pch = 19, cex = 0.3,
##   col = "red", add = TRUE)

## ----Torgegram, include=FALSE,cache=FALSE--------------------------------
mf04.Torg <- Torgegram(mf04p, "Summer_mn", nlag = 20, maxlag = 50000)
plot(mf04.Torg)

## ----eval=FALSE----------------------------------------------------------
## mf04.Torg <- Torgegram(mf04p, "Summer_mn", nlag = 20, maxlag = 50000)
## plot(mf04.Torg)

## ----GauModel0-----------------------------------------------------------
mf04.glmssn0 <- glmssn(Summer_mn ~ ELEV_DEM + SLOPE, mf04p,
   CorModels = NULL, use.nugget = TRUE)
summary(mf04.glmssn0)

## ----lm0-----------------------------------------------------------------
summary(lm(Summer_mn ~ ELEV_DEM + SLOPE, getSSNdata.frame(mf04p)))

## ----GauModel1, cache = TRUE---------------------------------------------
mf04.glmssn1 <- glmssn(Summer_mn ~ ELEV_DEM + SLOPE, mf04p,
   CorModels = c("Exponential.tailup", "Exponential.taildown",
      "Exponential.Euclid"), addfunccol = "afvArea")
summary(mf04.glmssn1)

## ----BinModel1, cache = TRUE---------------------------------------------
mf04.glmssnBin <- glmssn(MaxOver20 ~ ELEV_DEM + SLOPE, mf04p,
  CorModels = c("Mariah.tailup", "Spherical.taildown"),
  family = "binomial", addfunccol = "afvArea")
summary(mf04.glmssnBin)

## ----PoiModel1, cache = TRUE---------------------------------------------
mf04.glmssnPoi <- glmssn(C16 ~ ELEV_DEM + SLOPE, mf04p,
  CorModels = c("LinearSill.tailup", "LinearSill.taildown"),
  family = "poisson", addfunccol = "afvArea")
summary(mf04.glmssnPoi)

## ----Model1, include=FALSE, cache=FALSE----------------------------------
mf04.resid1 <- residuals(mf04.glmssn1)
names( getSSNdata.frame(mf04.resid1) )
plot(mf04.resid1)

## ------------------------------------------------------------------------
mf04.resid1 <- residuals(mf04.glmssn1)
names( getSSNdata.frame(mf04.resid1) )

## ----eval=FALSE----------------------------------------------------------
## plot(mf04.resid1)

## ----ResidHist,include=FALSE,cache=FALSE---------------------------------
par(mfrow = c(1, 2))
hist(mf04.resid1)
hist(mf04p, "Summer_mn")

## ----eval=FALSE----------------------------------------------------------
## par(mfrow = c(1, 2))
## hist(mf04.resid1)
## hist(mf04p, "Summer_mn")

## ------------------------------------------------------------------------
ObsDFr <- getSSNdata.frame(mf04.resid1)
ObsDF <- getSSNdata.frame(mf04p)
indOutlier <- ObsDFr["_resid_"] < -3
ObsDF[indOutlier, "Summer_mn"] <- NA
mf04c <- putSSNdata.frame(ObsDF, mf04p)

## ----cache=FALSE---------------------------------------------------------
mf04c.glmssn0 <- glmssn(Summer_mn ~ ELEV_DEM + SLOPE, mf04c,
   CorModels = c("Exponential.tailup", "Exponential.taildown",
   "Exponential.Euclid"), addfunccol = "afvArea", EstMeth = "ML")
summary(mf04c.glmssn0)

## ----cache=FALSE---------------------------------------------------------
mf04c.glmssn1 <- glmssn(Summer_mn ~ ELEV_DEM, mf04c,
   CorModels = c("Exponential.tailup", "Exponential.taildown"),
   addfunccol = "afvArea", EstMeth = "ML")
summary(mf04c.glmssn1)

## ----LOOCV, include=FALSE------------------------------------------------
cv.out <- CrossValidationSSN(mf04c.glmssn1)
par(mfrow = c(1, 2))
plot(mf04c.glmssn1$sampinfo$z,
   cv.out[, "cv.pred"], pch = 19,
   xlab = "Observed Data", ylab = "LOOCV Prediction")
abline(0, 1)
plot( na.omit( getSSNdata.frame(mf04c)[, "Summer_mn"]),
   cv.out[, "cv.se"], pch = 19,
   xlab = "Observed Data", ylab = "LOOCV Prediction SE")

## ----eval=FALSE----------------------------------------------------------
## cv.out <- CrossValidationSSN(mf04c.glmssn1)
## par(mfrow = c(1, 2))
## plot(mf04c.glmssn1$sampinfo$z,
##    cv.out[, "cv.pred"], pch = 19,
##    xlab = "Observed Data", ylab = "LOOCV Prediction")
## abline(0, 1)
## plot( na.omit( getSSNdata.frame(mf04c)[, "Summer_mn"]),
##    cv.out[, "cv.se"], pch = 19,
##    xlab = "Observed Data", ylab = "LOOCV Prediction SE")

## ----LOOCVSummary--------------------------------------------------------
CrossValidationStatsSSN(mf04c.glmssn1)

## ------------------------------------------------------------------------
GR2(mf04c.glmssn1)
varcomp(mf04c.glmssn1)

## ------------------------------------------------------------------------
AIC(mf04c.glmssn0)
AIC(mf04c.glmssn1)

## ------------------------------------------------------------------------
mf04c.glmssn1 <- glmssn(Summer_mn ~ ELEV_DEM, mf04c,
   CorModels = c("Exponential.tailup", "Exponential.taildown"),
   addfunccol = "afvArea")
mf04c.glmssn2 <- glmssn(Summer_mn ~ ELEV_DEM,  mf04c,
   CorModels = c("LinearSill.tailup", "Mariah.taildown"),
   addfunccol = "afvArea")
mf04c.glmssn3 <- glmssn(Summer_mn ~ ELEV_DEM , mf04c,
   CorModels =  c("Mariah.tailup", "LinearSill.taildown"),
   addfunccol = "afvArea")
mf04c.glmssn4 <- glmssn(Summer_mn ~ ELEV_DEM, mf04c,
   CorModels = c("Spherical.tailup", "Spherical.taildown"),
   addfunccol = "afvArea")
mf04c.glmssn5 <- glmssn(Summer_mn ~ ELEV_DEM, mf04c,
   CorModels = "Exponential.Euclid",
   addfunccol = "afvArea")

## ------------------------------------------------------------------------
options(digits = 4)
InfoCritCompare(list(mf04c.glmssn1, mf04c.glmssn2,
   mf04c.glmssn3, mf04c.glmssn4, mf04c.glmssn5))
options(digits = 7)

## ------------------------------------------------------------------------
summary(mf04c.glmssn2)

## ----Residuals,include=FALSE,cache=FALSE---------------------------------
mf04c.resid2 <- residuals(mf04c.glmssn2,
    cross.validation = TRUE)
mf04c.resid2.cv.std <-
    getSSNdata.frame(mf04c.resid2)[, "_resid.crossv_"] /
    getSSNdata.frame(mf04c.resid2)[, "_CrossValStdErr_"]
hist(mf04c.resid2.cv.std)

## ----eval=FALSE----------------------------------------------------------
## mf04c.resid2 <- residuals(mf04c.glmssn2,
##     cross.validation = TRUE)
## mf04c.resid2.cv.std <-
##     getSSNdata.frame(mf04c.resid2)[, "_resid.crossv_"] /
##     getSSNdata.frame(mf04c.resid2)[, "_CrossValStdErr_"]
## hist(mf04c.resid2.cv.std)

## ----TorgRes,include=FALSE,cache=FALSE-----------------------------------
plot(Torgegram(mf04c.resid2, "_resid_", nlag = 8, maxlag = 25000))

## ----eval=FALSE----------------------------------------------------------
## plot(Torgegram(mf04c.resid2, "_resid_", nlag = 8, maxlag = 25000))

## ----Preds1,include=FALSE,cache=FALSE------------------------------------
mf04c.pred1km <- predict(mf04c.glmssn4, "pred1km")
plot(mf04c.pred1km, SEcex.max = 1, SEcex.min = .5/3*2,
     breaktype = "user", brks = brks)

## ----eval=FALSE----------------------------------------------------------
## mf04c.pred1km <- predict(mf04c.glmssn4, "pred1km")
## plot(mf04c.pred1km, SEcex.max = 1, SEcex.min = .5/3*2,
##      breaktype = "user", brks = brks)

## ----Preds2,include=FALSE------------------------------------------------
plot(mf04c, "Summer_mn", pch = 1, cex = 3,
   xlab = "x-coordinate", ylab = "y-coordinate",
   xlim = c(-1511000,-1500000), ylim = c(2525000,2535000))
mf04c.glmssn4.Knapp <- predict(mf04c.glmssn4, "Knapp")
plot(mf04c.glmssn4.Knapp, "Summer_mn", add = TRUE,
   xlim = c(-1511000,-1500000), ylim = c(2525000,2535000))

## ----eval=FALSE----------------------------------------------------------
## plot(mf04c, "Summer_mn", pch = 1, cex = 3,
##    xlab = "x-coordinate", ylab = "y-coordinate",
##    xlim = c(-1511000,-1500000), ylim = c(2525000,2535000))
## mf04c.glmssn4.Knapp <- predict(mf04c.glmssn4, "Knapp")
## plot(mf04c.glmssn4.Knapp, "Summer_mn", add = TRUE,
##    xlim = c(-1511000,-1500000), ylim = c(2525000,2535000))

## ------------------------------------------------------------------------
mf04c.glmssn4.BPKnapp <- BlockPredict(mf04c.glmssn4, "Knapp")
mf04c.glmssn4.BPKnapp

## ------------------------------------------------------------------------
mf04c.glmssn4.BPCapeHorn <- BlockPredict(mf04c.glmssn4, "CapeHorn")
mf04c.glmssn4.BPCapeHorn

## ------------------------------------------------------------------------
mf04c.missingobs <- predict(mf04c.glmssn4, "_MissingObs_")
getPreds(mf04c.missingobs, pred.type = "pred")
with(getSSNdata.frame(mf04p), Summer_mn[pid==29])

## ----eval=FALSE----------------------------------------------------------
## createSSN(n, obsDesign, predDesign = noPoints, path,
##    importToR = FALSE, treeFunction = igraphKamadaKawai)

## ----SimIterative, include=FALSE, cache=FALSE----------------------------
set.seed(12)
iterative.ssn <- createSSN(n = c(30, 10),
   obsDesign = binomialDesign(c(10,10)),
   importToR = TRUE, path = "./SimIterative.ssn",
   treeFunction = iterativeTreeLayout)
plot(iterative.ssn, lwdLineCol = "addfunccol", lwdLineEx = 8,
   lineCol = "blue", cex = 2, xlab = "x-coordinate",
   ylab = "y-coordinate", pch = 1)

## ----eval=FALSE, cache=FALSE---------------------------------------------
## set.seed(12)
## iterative.ssn <- createSSN(n = c(30, 10),
##    obsDesign = binomialDesign(c(10,10)),
##    importToR = TRUE, path = "./SimIterative.ssn",
##    treeFunction = iterativeTreeLayout)
## plot(iterative.ssn, lwdLineCol = "addfunccol", lwdLineEx = 8,
##    lineCol = "blue", cex = 2, xlab = "x-coordinate",
##    ylab = "y-coordinate", pch = 1)

## ----SimSSN1, include=FALSE, cache=FALSE---------------------------------
set.seed(101)
raw.ssn <- createSSN(n = c(10, 10, 10),
   obsDesign = binomialDesign(c(40, 40, 40)),
   predDesign = systematicDesign(c(0.2, 0.4, 0.8)), importToR = TRUE,
   path = "./raw.ssn")
plot(raw.ssn, lwdLineCol = "addfunccol", lwdLineEx = 8,
   lineCol = "blue", cex = 2, xlab = "x-coordinate",
   ylab = "y-coordinate", pch = 1)
plot(raw.ssn, PredPointsID = "preds", add = TRUE, cex = .5, pch = 19,
   col = "green")

## ----eval=FALSE, cache=FALSE---------------------------------------------
## set.seed(101)
## raw.ssn <- createSSN(n = c(10, 10, 10),
##    obsDesign = binomialDesign(c(40, 40, 40)),
##    predDesign = systematicDesign(c(0.2, 0.4, 0.8)), importToR = TRUE,
##    path = "./raw.ssn")
## plot(raw.ssn, lwdLineCol = "addfunccol", lwdLineEx = 8,
##    lineCol = "blue", cex = 2, xlab = "x-coordinate",
##    ylab = "y-coordinate", pch = 1)
## plot(raw.ssn, PredPointsID = "preds", add = TRUE, cex = .5, pch = 19,
##    col = "green")

## ----SimHardcore, include=FALSE, cache=FALSE-----------------------------
set.seed(13)
hardcore.ssn <- createSSN(n = c(10, 10),
   obsDesign = hardCoreDesign(c(200, 200), c(0.2, 0.4)),
   importToR = TRUE, path = "./SimHardcore.ssn")
plot(hardcore.ssn, lwdLineCol = "addfunccol", lwdLineEx = 8,
   lineCol = "blue", cex = 2, xlab = "x-coordinate",
   ylab = "y-coordinate", pch = 1)
plot(hardcore.ssn, PredPointsID = NULL, add = TRUE, cex = .5,
   pch = 19, col = "green")

## ----eval=FALSE----------------------------------------------------------
## set.seed(13)
## hardcore.ssn <- createSSN(n = c(10, 10),
##    obsDesign = hardCoreDesign(c(200, 200), c(0.2, 0.4)),
##    importToR = TRUE, path = "./SimHardcore.ssn")
## plot(hardcore.ssn, lwdLineCol = "addfunccol", lwdLineEx = 8,
##    lineCol = "blue", cex = 2, xlab = "x-coordinate",
##    ylab = "y-coordinate", pch = 1)
## plot(hardcore.ssn, PredPointsID = NULL, add = TRUE, cex = .5,
##    pch = 19, col = "green")

## ----cache=FALSE---------------------------------------------------------
createDistMat(raw.ssn, "preds", o.write=TRUE, amongpred = TRUE)

## ------------------------------------------------------------------------
rawDFobs <- getSSNdata.frame(raw.ssn, "Obs")
rawDFpred <- getSSNdata.frame(raw.ssn, "preds")

## ------------------------------------------------------------------------
rawDFobs[,"X1"] <- rnorm(length(rawDFobs[,1]))
rawDFpred[,"X1"] <- rnorm(length(rawDFpred[,1]))
rawDFobs[,"X2"] <- rnorm(length(rawDFobs[,1]))
rawDFpred[,"X2"] <- rnorm(length(rawDFpred[,1]))

## ------------------------------------------------------------------------
rawDFobs[,"F1"] <- as.factor(sample.int(4,length(rawDFobs[,1]),
   replace = TRUE))
rawDFpred[,"F1"] <- as.factor(sample.int(4,length(rawDFpred[,1]),
   replace = TRUE))

## ------------------------------------------------------------------------
rawDFobs[,"RE1"] <- as.factor(sample(1:3,length(rawDFobs[,1]),
   replace = TRUE))
rawDFobs[,"RE2"] <- as.factor(sample(1:4,length(rawDFobs[,1]),
   replace = TRUE))
rawDFpred[,"RE1"] <- as.factor(sample(1:3,length(rawDFpred[,1]),
   replace = TRUE))
rawDFpred[,"RE2"] <- as.factor(sample(1:4,length(rawDFpred[,1]),
   replace = TRUE))

## ------------------------------------------------------------------------
names(rawDFobs)
names(rawDFpred)

## ----cache=FALSE---------------------------------------------------------
set.seed(102)
sim.out <- SimulateOnSSN(raw.ssn, ObsSimDF = rawDFobs,
   PredSimDF = rawDFpred, PredID = "preds",
   formula = ~ X1 + X2 + F1, coefficients = c(10,1,0,-2,0,2),
   CorModels = c("LinearSill.tailup", "Mariah.taildown",
   "Exponential.Euclid", "RE1", "RE2"), use.nugget = TRUE,
   CorParms = c(3, 10, 2, 10, 1, 5, 1, .5, .1),
   addfunccol = "addfunccol")

## ------------------------------------------------------------------------
with(rawDFobs, colnames(model.matrix( ~ X1 + X2 + F1)))

## ------------------------------------------------------------------------
sim.out$FixedEffects
sim.out$CorParms

## ------------------------------------------------------------------------
sim.ssn <- sim.out$ssn.object

## ----SimSSN2,include=FALSE,results='hide',cache=FALSE--------------------
plot(sim.ssn, "Sim_Values",
   xlab = "x-coordinate", ylab = "y-coordinate",
   cex = 1.5)

## ----eval=FALSE----------------------------------------------------------
## plot(sim.ssn, "Sim_Values",
##    xlab = "x-coordinate", ylab = "y-coordinate",
##    cex = 1.5)

## ------------------------------------------------------------------------
simDFobs <- getSSNdata.frame(sim.ssn, "Obs")
simDFpred <- getSSNdata.frame(sim.ssn, "preds")

## ------------------------------------------------------------------------
simpreds <- simDFpred[,"Sim_Values"]
simDFpred[,"Sim_Values"] <- NA
sim.ssn <- putSSNdata.frame(simDFpred, sim.ssn, "preds")

## ----cache=FALSE---------------------------------------------------------
glmssn.out <- glmssn(Sim_Values ~ X1 + X2 + F1, sim.ssn,
   CorModels = c("LinearSill.tailup", "Mariah.taildown",
   "Exponential.Euclid", "RE1", "RE2"),
   addfunccol = "addfunccol")

## ------------------------------------------------------------------------
summary(glmssn.out)

## ----SimTvP,include=FALSE,cache=FALSE------------------------------------
glmssn.pred <- predict(glmssn.out,"preds")
predDF <- getSSNdata.frame(glmssn.pred, "preds")
plot(simpreds, predDF[,"Sim_Values"], xlab = "True",
   ylab = "Predicted", pch = 19)

## ----eval=FALSE,cache=FALSE----------------------------------------------
## glmssn.pred <- predict(glmssn.out,"preds")
## predDF <- getSSNdata.frame(glmssn.pred, "preds")
## plot(simpreds, predDF[,"Sim_Values"], xlab = "True",
##    ylab = "Predicted", pch = 19)

