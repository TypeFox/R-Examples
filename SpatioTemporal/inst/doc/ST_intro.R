### R code from vignette source 'ST_intro.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: ST_intro.Rnw:124-128
###################################################
##initial set-up: load libraries
library(SpatioTemporal)
library(plotrix) 
library(maps)


###################################################
### code chunk number 2: ST_intro.Rnw:678-723
###################################################
##load precomputed results
data(mesa.data.raw, package="SpatioTemporal")
data(mesa.model, package="SpatioTemporal")
data(est.mesa.model, package="SpatioTemporal")
  
##create model without spatiotemporal covariate 
##(needed to be able to add time-points)
mesa.data <- createSTdata(mesa.data.raw$obs, mesa.data.raw$X, n.basis=2,
                          SpatioTemporal=mesa.data.raw["lax.conc.1500"])
mesa.model <- createSTmodel(mesa.data, LUR=mesa.model$LUR.list, 
                            ST=NULL, cov.beta=mesa.model$cov.beta, 
                            cov.nu=mesa.model$cov.nu,
                            locations=mesa.model$locations.list)

##restrict mesa.data to one location
mesa.model.obs <- dropObservations(mesa.model,
                                   mesa.model$obs$ID == "60590007")
  
##restrict mesa.data to one location
mesa.data$covars <- mesa.data$covars[mesa.data$covars$ID=="60590007",,drop=FALSE]
mesa.data$obs <- mesa.data$obs[mesa.data$obs$ID=="60590007",,drop=FALSE]
mesa.data$SpatioTemporal <- NULL
##add more temporal points (every 3.5 days)
mesa.data.2 <- updateTrend(mesa.data, fnc=mesa.data$trend.fnc,
                           extra.dates=seq(min(mesa.data$trend$date),
                             max(mesa.data$trend$date), by=3.5))
  
##nugget (using first value is fine; it is an AQS-site)
pars <- coef(est.mesa.model, pars="cov")$par
nugget <- loglikeSTgetPars(pars, mesa.model)$cov.nu$nugget[1]
  
##predictions at this location (every 14 days)
pred <- predict(mesa.model.obs, pars, mesa.data,
                nugget.unobs=nugget)
##new predictions
pred.2 <- predict(mesa.model.obs, pars,  mesa.data.2,
                  nugget.unobs=nugget)

##predictions for log-field
pred.log <- predict(mesa.model.obs, pars, mesa.data,
                    nugget.unobs=nugget, transform="unbiased", 
                    type="p")
pred.log.2 <- predict(mesa.model.obs, pars, mesa.data.2,
                      nugget.unobs=nugget, transform="unbiased", 
                      type="p")


###################################################
### code chunk number 3: figPredictExample (eval = FALSE)
###################################################
## ##prediction interval for additional time-points
## plot(pred.2, ID="60590007", main="Predictions at AQS 60590007",
##      xlab="", ylab="NOx (log pbb)", xaxt="n",
##      xlim=as.Date(c("2006-07-01","2007-01-01")),
##      ylim=c(2.3,5.2), pch=NA, pred.var=TRUE,
##      lty=NA, col=c(1,1,"lightgrey"))
## ##prediction interval for original timepoints
## plot(pred, ID="60590007", pred.var=TRUE, add=TRUE, 
##      lty=NA, pch=NA, col=c(1,1,"darkgrey"))
## ##confidence interval for original timepoints
## plot(pred, ID="60590007", add=TRUE, lty=NA, pch=NA, col=c(1,1,"white"))
## ##predictions for additional time-points, due to mean component
## plot(pred.2, ID="60590007", add=TRUE, pred.type = "EX.mu", 
##      lty=c(4,NA), lwd=2, pch=NA, col=c(1,NA,NA))
## ##predictions for additional time-points
## plot(pred.2, ID="60590007", add=TRUE, pred.type = "EX", 
##      lty=c(1,NA), lwd=2, pch=NA, col=c(1,NA,NA))
## ##predictions for additional time-points, due to beta-field
## plot(pred.2, ID="60590007", add=TRUE, pred.type = "EX.mu.beta", 
##      lty=c(2,NA), lwd=2, pch=NA, col=c("green",NA,NA))
## ##observations and predictions at original time-scale
## plot(pred, ID="60590007", STmodel=mesa.data.2, add=TRUE, 
##      pred.type = "EX", lty=NA, pch=c(3,4), col=c(1,2,NA))


###################################################
### code chunk number 4: figPredictLogExample (eval = FALSE)
###################################################
## ##prediction interval for additional time-points
## plot(pred.log.2, pred.type="EX.pred", ID="60590007", 
##      main="",  xlab="", ylab="NOx (pbb)", 
##      xlim=as.Date(c("2006-07-01","2007-01-01")),
##      ylim=c(0,185), pred.var=TRUE, pch=NA, 
##      lty=NA, col=c(1,1,"lightgrey"))
## ##prediction interval for original timepoints
## plot(pred.log, pred.type="EX.pred", ID="60590007", pred.var=TRUE, add=TRUE,
##      lty=NA, pch=NA, col=c(1,1,"darkgrey"))
## ##confidence interval for original timepoints
## plot(pred.log, ID="60590007", add=TRUE, 
##      lty=NA, pch=NA, col=c(1,1,"white"))
## ##predictions for additional time-points, due to mean component
## plot(pred.log.2, ID="60590007", add=TRUE, pred.type = "EX.mu", 
##      lty=c(4,NA), lwd=2, pch=NA, col=c(1,NA,NA))
## ##predictions for additional time-points
## plot(pred.log.2, ID="60590007", add=TRUE, pred.type = "EX.pred", 
##      lty=c(1,NA), pch=NA, col=c("orange",NA,NA), lwd=3, pred.var=TRUE)
## plot(pred.log.2, ID="60590007", add=TRUE, pred.type = "EX", 
##      lty=c(1,NA), pch=NA, col=c(1,NA,NA), lwd=2)
## ##predictions for additional time-points, due to beta-field
## plot(pred.log.2, ID="60590007", add=TRUE, pred.type = "EX.mu.beta", 
##      lty=c(2,NA), lwd=2, pch=NA, col=c("green",NA,NA))
## ##observations and predictions at original time-scale
## plot(pred.log, ID="60590007", STmodel=mesa.data.2, add=TRUE, 
##      pred.type = "EX", lty=NA, pch=c(3,4), col=c(1,2,NA))
## legend("topleft", c("Observations", 
##                     "Predictions",
##                     "Predictions, incl. nugget",
##                     "Contribution from beta",
##                     "Contribution from mean",
##                     "95% CI"), bty="n",
##        lty=c(NA,1,1,2,4,NA), lwd=c(NA,2,2,2,2,NA),
##        pch=c(4,3,NA,NA,NA,15), pt.cex=c(1,1,NA,NA,NA,2.5),
##        col=c("red", 1, "orange", "green", 1, "grey"))


###################################################
### code chunk number 5: ST_intro.Rnw:830-834
###################################################
par(mfcol=c(2,1), mar=c(.5,3.3,2,1), mgp=c(2,1,0))
##prediction interval for additional time-points
plot(pred.2, ID="60590007", main="Predictions at AQS 60590007",
     xlab="", ylab="NOx (log pbb)", xaxt="n",
     xlim=as.Date(c("2006-07-01","2007-01-01")),
     ylim=c(2.3,5.2), pch=NA, pred.var=TRUE,
     lty=NA, col=c(1,1,"lightgrey"))
##prediction interval for original timepoints
plot(pred, ID="60590007", pred.var=TRUE, add=TRUE, 
     lty=NA, pch=NA, col=c(1,1,"darkgrey"))
##confidence interval for original timepoints
plot(pred, ID="60590007", add=TRUE, lty=NA, pch=NA, col=c(1,1,"white"))
##predictions for additional time-points, due to mean component
plot(pred.2, ID="60590007", add=TRUE, pred.type = "EX.mu", 
     lty=c(4,NA), lwd=2, pch=NA, col=c(1,NA,NA))
##predictions for additional time-points
plot(pred.2, ID="60590007", add=TRUE, pred.type = "EX", 
     lty=c(1,NA), lwd=2, pch=NA, col=c(1,NA,NA))
##predictions for additional time-points, due to beta-field
plot(pred.2, ID="60590007", add=TRUE, pred.type = "EX.mu.beta", 
     lty=c(2,NA), lwd=2, pch=NA, col=c("green",NA,NA))
##observations and predictions at original time-scale
plot(pred, ID="60590007", STmodel=mesa.data.2, add=TRUE, 
     pred.type = "EX", lty=NA, pch=c(3,4), col=c(1,2,NA))
par(mar=c(2.3,3.3,0,1))
##prediction interval for additional time-points
plot(pred.log.2, pred.type="EX.pred", ID="60590007", 
     main="",  xlab="", ylab="NOx (pbb)", 
     xlim=as.Date(c("2006-07-01","2007-01-01")),
     ylim=c(0,185), pred.var=TRUE, pch=NA, 
     lty=NA, col=c(1,1,"lightgrey"))
##prediction interval for original timepoints
plot(pred.log, pred.type="EX.pred", ID="60590007", pred.var=TRUE, add=TRUE,
     lty=NA, pch=NA, col=c(1,1,"darkgrey"))
##confidence interval for original timepoints
plot(pred.log, ID="60590007", add=TRUE, 
     lty=NA, pch=NA, col=c(1,1,"white"))
##predictions for additional time-points, due to mean component
plot(pred.log.2, ID="60590007", add=TRUE, pred.type = "EX.mu", 
     lty=c(4,NA), lwd=2, pch=NA, col=c(1,NA,NA))
##predictions for additional time-points
plot(pred.log.2, ID="60590007", add=TRUE, pred.type = "EX.pred", 
     lty=c(1,NA), pch=NA, col=c("orange",NA,NA), lwd=3, pred.var=TRUE)
plot(pred.log.2, ID="60590007", add=TRUE, pred.type = "EX", 
     lty=c(1,NA), pch=NA, col=c(1,NA,NA), lwd=2)
##predictions for additional time-points, due to beta-field
plot(pred.log.2, ID="60590007", add=TRUE, pred.type = "EX.mu.beta", 
     lty=c(2,NA), lwd=2, pch=NA, col=c("green",NA,NA))
##observations and predictions at original time-scale
plot(pred.log, ID="60590007", STmodel=mesa.data.2, add=TRUE, 
     pred.type = "EX", lty=NA, pch=c(3,4), col=c(1,2,NA))
legend("topleft", c("Observations", 
                    "Predictions",
                    "Predictions, incl. nugget",
                    "Contribution from beta",
                    "Contribution from mean",
                    "95% CI"), bty="n",
       lty=c(NA,1,1,2,4,NA), lwd=c(NA,2,2,2,2,NA),
       pch=c(4,3,NA,NA,NA,15), pt.cex=c(1,1,NA,NA,NA,2.5),
       col=c("red", 1, "orange", "green", 1, "grey"))


###################################################
### code chunk number 6: ST_intro.Rnw:1176-1180
###################################################
library(SpatioTemporal)
library(plotrix) 
library(maps)
data(mesa.data.raw, package="SpatioTemporal")


###################################################
### code chunk number 7: ST_intro.Rnw:1188-1190
###################################################
mesa.data <- createSTdata(obs=mesa.data.raw$obs, covars=mesa.data.raw$X, 
   SpatioTemporal=list(lax.conc.1500=mesa.data.raw$lax.conc.1500))


###################################################
### code chunk number 8: ST_intro.Rnw:1208-1209
###################################################
table(mesa.data.raw$X$type)


###################################################
### code chunk number 9: figPlotSTdata (eval = FALSE)
###################################################
## layout(matrix(c(1,2,1,3), 2, 2))
## par(mar=c(2.3,3.3,2,1), mgp=c(2,1,0))
## plot(mesa.data, "loc", main="Occurrence of Observations", xlab="", 
##      ylab="Location", col=c("black", "red"), legend.loc=NULL)
## par(mar=c(3.3,3.3,2,1))
## qqnorm(mesa.data, line=1)
## scatterPlot(mesa.data, covar="km.to.coast", xlab="Distance to coast",
##             ylab="NOx (log ppb)", pch=19, cex=.25,
##             smooth.args=list(span=4/5,degree=2))


###################################################
### code chunk number 10: ST_intro.Rnw:1236-1237
###################################################
layout(matrix(c(1,2,1,3), 2, 2))
par(mar=c(2.3,3.3,2,1), mgp=c(2,1,0))
plot(mesa.data, "loc", main="Occurrence of Observations", xlab="", 
     ylab="Location", col=c("black", "red"), legend.loc=NULL)
par(mar=c(3.3,3.3,2,1))
qqnorm(mesa.data, line=1)
scatterPlot(mesa.data, covar="km.to.coast", xlab="Distance to coast",
            ylab="NOx (log ppb)", pch=19, cex=.25,
            smooth.args=list(span=4/5,degree=2))


###################################################
### code chunk number 11: ST_intro.Rnw:1256-1257
###################################################
D <- createDataMatrix(mesa.data)


###################################################
### code chunk number 12: ST_intro.Rnw:1268-1269
###################################################
SVD.cv <- SVDsmoothCV(D, 0:4)


###################################################
### code chunk number 13: ST_intro.Rnw:1271-1272 (eval = FALSE)
###################################################
## SVD.cv <- SVDsmoothCV(D, 0:4)


###################################################
### code chunk number 14: ST_intro.Rnw:1277-1278
###################################################
print(SVD.cv)


###################################################
### code chunk number 15: figSVDcv (eval = FALSE)
###################################################
## plot(SVD.cv)


###################################################
### code chunk number 16: ST_intro.Rnw:1294-1296
###################################################
par(mgp=c(2,1,0))
plot(SVD.cv)


###################################################
### code chunk number 17: ST_intro.Rnw:1304-1305
###################################################
mesa.data <- updateTrend(mesa.data, n.basis=2)


###################################################
### code chunk number 18: ST_intro.Rnw:1310-1311
###################################################
smooth.trend <- calcSmoothTrends(mesa.data, n.basis=2, cv=TRUE)


###################################################
### code chunk number 19: ST_intro.Rnw:1313-1314 (eval = FALSE)
###################################################
## smooth.trend <- calcSmoothTrends(mesa.data, n.basis=2, cv=TRUE)


###################################################
### code chunk number 20: figsmooth_trends_CV (eval = FALSE)
###################################################
## mesa.data.cv <- vector("list",  length(smooth.trend$trend.fnc.cv))
## for(i in 1:length(mesa.data.cv)){
##   suppressMessages(mesa.data.cv[[i]] <- updateTrend(mesa.data, 
##                    fnc=smooth.trend$trend.fnc.cv[[i]]))
## }
## plot(mesa.data, main="Possible temporal trends",
##      xlab="", ylab="NOx (log ppb)", pch=c(19,NA), cex=.25)
## for(i in 1:length(mesa.data.cv)){
##   plot(mesa.data.cv[[i]], add=TRUE, col=i, pch=NA, lty=c(NA,2))
## }


###################################################
### code chunk number 21: ST_intro.Rnw:1335-1337
###################################################
par(mar=c(2.1,3.3,2,1), mgp=c(2,1,0))
mesa.data.cv <- vector("list",  length(smooth.trend$trend.fnc.cv))
for(i in 1:length(mesa.data.cv)){
  suppressMessages(mesa.data.cv[[i]] <- updateTrend(mesa.data, 
                   fnc=smooth.trend$trend.fnc.cv[[i]]))
}
plot(mesa.data, main="Possible temporal trends",
     xlab="", ylab="NOx (log ppb)", pch=c(19,NA), cex=.25)
for(i in 1:length(mesa.data.cv)){
  plot(mesa.data.cv[[i]], add=TRUE, col=i, pch=NA, lty=c(NA,2))
}


###################################################
### code chunk number 22: figsmooth_trends (eval = FALSE)
###################################################
## par(mar=c(3.3,3.3,1.5,1), mgp=c(2,1,0))
## layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow=TRUE))
## plot(mesa.data, "obs", ID="60370113", 
##      xlab="", ylab="NOx (log ppb)",
##      main="Temporal trend 60370113")
## plot(mesa.data, "res", ID="60370113", 
##      xlab="", ylab="NOx (log ppb)")
## plot(mesa.data, "acf", ID="60370113")
## plot(mesa.data, "pacf", ID="60370113")


###################################################
### code chunk number 23: ST_intro.Rnw:1369-1370
###################################################
par(mar=c(3.3,3.3,1.5,1), mgp=c(2,1,0))
layout(matrix(c(1,1,2,2,3,4), 3, 2, byrow=TRUE))
plot(mesa.data, "obs", ID="60370113", 
     xlab="", ylab="NOx (log ppb)",
     main="Temporal trend 60370113")
plot(mesa.data, "res", ID="60370113", 
     xlab="", ylab="NOx (log ppb)")
plot(mesa.data, "acf", ID="60370113")
plot(mesa.data, "pacf", ID="60370113")


###################################################
### code chunk number 24: figDeterministicTrends (eval = FALSE)
###################################################
## mesa.data.fnc <- updateTrend(mesa.data, fnc=function(x){
##   x = 2*pi*as.numeric(x)/365; 
##   return( cbind(x, sin(x), cos(x)) )})
## par(mfrow=c(2,1), mar=c(2.3,3.3,1.5,1), mgp=c(2,1,0))
## for(i in c("60370016","60371103")){
##   plot(mesa.data, ID=i, pch=c(19,NA), cex=.25, xlab="",
##        ylab="NOx (log ppb)", main=paste("AQS site",i))
##   plot(mesa.data.fnc, ID=i, add=TRUE, col=2, pch=NA)
## }


###################################################
### code chunk number 25: ST_intro.Rnw:1400-1401
###################################################
mesa.data.fnc <- updateTrend(mesa.data, fnc=function(x){
  x = 2*pi*as.numeric(x)/365; 
  return( cbind(x, sin(x), cos(x)) )})
par(mfrow=c(2,1), mar=c(2.3,3.3,1.5,1), mgp=c(2,1,0))
for(i in c("60370016","60371103")){
  plot(mesa.data, ID=i, pch=c(19,NA), cex=.25, xlab="",
       ylab="NOx (log ppb)", main=paste("AQS site",i))
  plot(mesa.data.fnc, ID=i, add=TRUE, col=2, pch=NA)
}


###################################################
### code chunk number 26: figBetaLm (eval = FALSE)
###################################################
## beta.lm <- estimateBetaFields(mesa.data)
## par(mfrow=c(1,2), mar=c(3.3,2.3,1.5,1), mgp=c(2,1,0))
## plotCI(mesa.data$covars$log10.m.to.a1, beta.lm$beta[,1], 
##        uiw=1.96*beta.lm$beta.sd[,1], ylab="", xlab="Distance to A1-road",
##        main="Beta-field for f1(t)")
## plotCI(mesa.data$covars$km.to.coast, beta.lm$beta[,2], 
##        uiw=1.96*beta.lm$beta.sd[,2], ylab="", xlab="Distance to coast",
##        main="Beta-field for f2(t)")


###################################################
### code chunk number 27: ST_intro.Rnw:1443-1444
###################################################
beta.lm <- estimateBetaFields(mesa.data)
par(mfrow=c(1,2), mar=c(3.3,2.3,1.5,1), mgp=c(2,1,0))
plotCI(mesa.data$covars$log10.m.to.a1, beta.lm$beta[,1], 
       uiw=1.96*beta.lm$beta.sd[,1], ylab="", xlab="Distance to A1-road",
       main="Beta-field for f1(t)")
plotCI(mesa.data$covars$km.to.coast, beta.lm$beta[,2], 
       uiw=1.96*beta.lm$beta.sd[,2], ylab="", xlab="Distance to coast",
       main="Beta-field for f2(t)")


###################################################
### code chunk number 28: ST_intro.Rnw:1458-1461
###################################################
LUR <- list(~log10.m.to.a1+s2000.pop.div.10000+km.to.coast,
            ~km.to.coast, ~km.to.coast)
cov.beta <- list(covf="exp", nugget=FALSE)


###################################################
### code chunk number 29: ST_intro.Rnw:1465-1466
###################################################
cov.nu <- list(covf="exp", nugget=~type, random.effect=FALSE)


###################################################
### code chunk number 30: ST_intro.Rnw:1495-1500
###################################################
locations <- list(coords=c("x","y"), long.lat=c("long","lat"), 
                  others="type")
mesa.model <- createSTmodel(mesa.data, LUR=LUR, ST="lax.conc.1500",
                            cov.beta=cov.beta, cov.nu=cov.nu,
                            locations=locations)


###################################################
### code chunk number 31: ST_intro.Rnw:1526-1530
###################################################
dim <- loglikeSTdim(mesa.model)
x.init <- cbind(c( rep(2, dim$nparam.cov-1), 0),
                c( rep(c(1,-3), dim$m+1), -3, 0))
rownames(x.init) <- loglikeSTnames(mesa.model, all=FALSE)


###################################################
### code chunk number 32: ST_intro.Rnw:1533-1534 (eval = FALSE)
###################################################
## est.mesa.model <- estimate(mesa.model, x.init, type="p", hessian.all=TRUE)


###################################################
### code chunk number 33: ST_intro.Rnw:1541-1543
###################################################
data(est.mesa.model, package="SpatioTemporal")
print(est.mesa.model)


###################################################
### code chunk number 34: ST_intro.Rnw:1567-1570
###################################################
pred <- predict(mesa.model, est.mesa.model, LTA=TRUE, type="p")
pred.log <- predict(mesa.model, est.mesa.model, LTA=TRUE,
                    transform="unbiased", type="p")


###################################################
### code chunk number 35: ST_intro.Rnw:1572-1575 (eval = FALSE)
###################################################
## pred <- predict(mesa.model, est.mesa.model, LTA=TRUE, type="p")
## pred.log <- predict(mesa.model, est.mesa.model, LTA=TRUE,
##                     transform="unbiased", type="p")


###################################################
### code chunk number 36: figBetaEX (eval = FALSE)
###################################################
## par(mfrow=c(2,2), mar=c(3.3,3.3,1.5,1), mgp=c(2,1,0), pty="s")
## for(i in 1:3){
##   plotCI(x=beta.lm$beta[,i], y=pred$beta$EX[,i], 
##          uiw=1.96*beta.lm$beta.sd[,i], err="x", 
##          pch=NA, sfrac=0.005,
##          main=paste("Beta-field for f", i, "(t)", sep=""),
##          xlab="Empirical estimate",
##          ylab="Spatio-Temporal Model")
##   plotCI(x=beta.lm$beta[,i], y=pred$beta$EX[,i], 
##          uiw=1.96*sqrt(pred$beta$VX[,i]),
##          add=TRUE, pch=NA, sfrac=0.005)
##   abline(0, 1, col="grey")
## }


###################################################
### code chunk number 37: ST_intro.Rnw:1611-1612
###################################################
par(mfrow=c(2,2), mar=c(3.3,3.3,1.5,1), mgp=c(2,1,0), pty="s")
for(i in 1:3){
  plotCI(x=beta.lm$beta[,i], y=pred$beta$EX[,i], 
         uiw=1.96*beta.lm$beta.sd[,i], err="x", 
         pch=NA, sfrac=0.005,
         main=paste("Beta-field for f", i, "(t)", sep=""),
         xlab="Empirical estimate",
         ylab="Spatio-Temporal Model")
  plotCI(x=beta.lm$beta[,i], y=pred$beta$EX[,i], 
         uiw=1.96*sqrt(pred$beta$VX[,i]),
         add=TRUE, pch=NA, sfrac=0.005)
  abline(0, 1, col="grey")
}


###################################################
### code chunk number 38: figLTA (eval = FALSE)
###################################################
## par(mfrow=c(1,2), mar=c(3.3,3.3,1.5,1), mgp=c(2,1,0))
## with(pred$LTA, plotCI(colMeans(D, na.rm=TRUE), EX, uiw=1.96*sqrt(VX.pred),
##                       xlim=c(2.9,4.4), ylim=c(2.9,4.4),
##                       ylab="Predictions", xlab="Observations",
##                       main="Average NOx (log ppb)"))
## abline(0, 1, col="grey")
## with(pred.log$LTA, plotCI(colMeans(exp(D), na.rm=TRUE), 
##                           EX, uiw=1.96*sqrt(VX.pred),
##                           xlim=c(25,95), ylim=c(25,95),
##                           ylab="Predictions", xlab="Observations",
##                           main="Average NOx (ppb)"))
## abline(0, 1, col="grey")


###################################################
### code chunk number 39: ST_intro.Rnw:1642-1643
###################################################
par(mfrow=c(1,2), mar=c(3.3,3.3,1.5,1), mgp=c(2,1,0))
with(pred$LTA, plotCI(colMeans(D, na.rm=TRUE), EX, uiw=1.96*sqrt(VX.pred),
                      xlim=c(2.9,4.4), ylim=c(2.9,4.4),
                      ylab="Predictions", xlab="Observations",
                      main="Average NOx (log ppb)"))
abline(0, 1, col="grey")
with(pred.log$LTA, plotCI(colMeans(exp(D), na.rm=TRUE), 
                          EX, uiw=1.96*sqrt(VX.pred),
                          xlim=c(25,95), ylim=c(25,95),
                          ylab="Predictions", xlab="Observations",
                          main="Average NOx (ppb)"))
abline(0, 1, col="grey")


###################################################
### code chunk number 40: ST_intro.Rnw:1658-1659
###################################################
Ind.cv <- createCV(mesa.model, groups=10, min.dist=.1)


###################################################
### code chunk number 41: ST_intro.Rnw:1666-1668
###################################################
ID.cv <- sapply(split(mesa.model$obs$ID, Ind.cv), unique)
print( sapply(ID.cv, length) )


###################################################
### code chunk number 42: ST_intro.Rnw:1671-1672
###################################################
table(Ind.cv)


###################################################
### code chunk number 43: ST_intro.Rnw:1676-1677
###################################################
print(mesa.model$D.beta[ID.cv[[10]],ID.cv[[10]]])


###################################################
### code chunk number 44: ST_intro.Rnw:1692-1693
###################################################
x.init <- coef(est.mesa.model, pars="cov")[,c("par","init")]


###################################################
### code chunk number 45: ST_intro.Rnw:1696-1697 (eval = FALSE)
###################################################
## est.cv.mesa <- estimateCV(mesa.model, x.init, Ind.cv)


###################################################
### code chunk number 46: ST_intro.Rnw:1700-1701
###################################################
data(est.cv.mesa, package="SpatioTemporal")


###################################################
### code chunk number 47: figParEstCV (eval = FALSE)
###################################################
## par(mfrow=c(1,1), mar=c(13.5,2.5,.5,.5), las=2)
## with(coef(est.mesa.model, pars="all"), 
##      plotCI((1:length(par))+.3, par, uiw=1.96*sd, 
##             col=2, xlab="", xaxt="n", ylab=""))
## boxplot(est.cv.mesa, "all", boxwex=.4, col="grey", add=TRUE)


###################################################
### code chunk number 48: ST_intro.Rnw:1716-1717
###################################################
par(mfrow=c(1,1), mar=c(13.5,2.5,.5,.5), las=2)
with(coef(est.mesa.model, pars="all"), 
     plotCI((1:length(par))+.3, par, uiw=1.96*sd, 
            col=2, xlab="", xaxt="n", ylab=""))
boxplot(est.cv.mesa, "all", boxwex=.4, col="grey", add=TRUE)


###################################################
### code chunk number 49: ST_intro.Rnw:1728-1732
###################################################
##pred.cv.mesa <- predictCV(mesa.model, est.cv.mesa, LTA=TRUE)
data(pred.cv.mesa, package="SpatioTemporal")
pred.cv.mesa.log <- predictCV(mesa.model, est.cv.mesa, 
                         LTA=TRUE, transform="unbiased")


###################################################
### code chunk number 50: ST_intro.Rnw:1745-1746
###################################################
summary(pred.cv.mesa.log)


###################################################
### code chunk number 51: figPredCV (eval = FALSE)
###################################################
## par(mar=c(3.3,3.3,1.5,1), mgp=c(2,1,0))
## layout(matrix(c(1,1,2,2,3,4,5,6), 4, 2, byrow=TRUE))
## plot(pred.cv.mesa, ID="60371601", xlab="", ylab="NOx (log ppb)", 
##      main="Predictions for 60371601", lty=c(1,NA), lwd=2, 
##      pch=c(NA,19), cex=.75)
## plot(pred.cv.mesa, ID="60371601", pred.type="EX.mu", 
##      lty=4, lwd=2, col="blue", add=TRUE)
## plot(pred.cv.mesa, ID="60371601", pred.type="EX.mu.beta", 
##      lty=2, lwd=2, col="green", add=TRUE)
## 
## plot(pred.cv.mesa.log, ID="60371601", xlab="", ylab="NOx (ppb)", 
##      main="Predictions for 60371601", pred.type="EX.pred", 
##      lty=c(1,NA), lwd=2, pch=c(NA,19), cex=.75)
## plot(pred.cv.mesa.log, ID="60371601", pred.type="EX.mu", 
##      lty=4, lwd=2, col="blue", add=TRUE)
## plot(pred.cv.mesa.log, ID="60371601", pred.type="EX.mu.beta", 
##      lty=2, lwd=2, col="green", add=TRUE)
## legend("topright", c("Observations", "Predictions",
##                      "Contribution from beta",
##                      "Contribution from mean",
##                      "95% CI"), bty="n",
##        lty=c(NA,1,2,4,NA), lwd=c(NA,2,2,2,NA),
##        pch=c(19,NA,NA,NA,15), pt.cex=c(.75,NA,NA,NA,2.5),
##        col=c("red", "black", "green", "blue", "grey"))
## 
## plot(pred.cv.mesa, "obs", ID="all", pch=c(19,NA), cex=.25, lty=c(NA,2), 
##      col=c("ID", "black", "grey"), xlab="Observations", 
##      ylab="Predictions", main="Cross-validation NOx (log ppb)")
## 
## with(pred.cv.mesa.log$pred.LTA, plotCI(obs, EX.pred, uiw=1.96*sqrt(VX.pred),
##                                   xlab="Observations", ylab="Predictions",
##                                   main="Temporal average NOx (ppb)"))
## abline(0, 1, col="grey")
## 
## I.season <- as.factor(as.POSIXlt(pred.cv.mesa$pred.obs$date)$mon+1)
## levels(I.season) <- c(rep("Winter",2), rep("Spring",3), 
##                       rep("Summer",3), rep("Fall",3), "Winter") 
## qqnorm(pred.cv.mesa, norm=TRUE, main="Normalised residuals",
##        col=I.season)
## legend("bottomright", legend=as.character(levels(I.season)),
##        pch=1, col=1:nlevels(I.season), bty="n")
## 
## scatterPlot(pred.cv.mesa, STdata=mesa.model, covar="log10.m.to.a1", 
##             group=I.season, col=c(2:5,1), type="res",
##             xlab="Distance to A1-Road", ylab="Residuals", 
##             main="Residuals (log ppb)")


###################################################
### code chunk number 52: ST_intro.Rnw:1844-1845
###################################################
par(mar=c(3.3,3.3,1.5,1), mgp=c(2,1,0))
layout(matrix(c(1,1,2,2,3,4,5,6), 4, 2, byrow=TRUE))
plot(pred.cv.mesa, ID="60371601", xlab="", ylab="NOx (log ppb)", 
     main="Predictions for 60371601", lty=c(1,NA), lwd=2, 
     pch=c(NA,19), cex=.75)
plot(pred.cv.mesa, ID="60371601", pred.type="EX.mu", 
     lty=4, lwd=2, col="blue", add=TRUE)
plot(pred.cv.mesa, ID="60371601", pred.type="EX.mu.beta", 
     lty=2, lwd=2, col="green", add=TRUE)

plot(pred.cv.mesa.log, ID="60371601", xlab="", ylab="NOx (ppb)", 
     main="Predictions for 60371601", pred.type="EX.pred", 
     lty=c(1,NA), lwd=2, pch=c(NA,19), cex=.75)
plot(pred.cv.mesa.log, ID="60371601", pred.type="EX.mu", 
     lty=4, lwd=2, col="blue", add=TRUE)
plot(pred.cv.mesa.log, ID="60371601", pred.type="EX.mu.beta", 
     lty=2, lwd=2, col="green", add=TRUE)
legend("topright", c("Observations", "Predictions",
                     "Contribution from beta",
                     "Contribution from mean",
                     "95% CI"), bty="n",
       lty=c(NA,1,2,4,NA), lwd=c(NA,2,2,2,NA),
       pch=c(19,NA,NA,NA,15), pt.cex=c(.75,NA,NA,NA,2.5),
       col=c("red", "black", "green", "blue", "grey"))

plot(pred.cv.mesa, "obs", ID="all", pch=c(19,NA), cex=.25, lty=c(NA,2), 
     col=c("ID", "black", "grey"), xlab="Observations", 
     ylab="Predictions", main="Cross-validation NOx (log ppb)")

with(pred.cv.mesa.log$pred.LTA, plotCI(obs, EX.pred, uiw=1.96*sqrt(VX.pred),
                                  xlab="Observations", ylab="Predictions",
                                  main="Temporal average NOx (ppb)"))
abline(0, 1, col="grey")

I.season <- as.factor(as.POSIXlt(pred.cv.mesa$pred.obs$date)$mon+1)
levels(I.season) <- c(rep("Winter",2), rep("Spring",3), 
                      rep("Summer",3), rep("Fall",3), "Winter") 
qqnorm(pred.cv.mesa, norm=TRUE, main="Normalised residuals",
       col=I.season)
legend("bottomright", legend=as.character(levels(I.season)),
       pch=1, col=1:nlevels(I.season), bty="n")

scatterPlot(pred.cv.mesa, STdata=mesa.model, covar="log10.m.to.a1", 
            group=I.season, col=c(2:5,1), type="res",
            xlab="Distance to A1-Road", ylab="Residuals", 
            main="Residuals (log ppb)")


