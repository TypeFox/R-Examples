### R code from vignette source 'ST_tutorial.Rnw'

###################################################
### code chunk number 1: ST_tutorial.Rnw:35-37
###################################################
library(SpatioTemporal)
data(mesa.data.raw, package="SpatioTemporal")


###################################################
### code chunk number 2: ST_tutorial.Rnw:53-54
###################################################
options(width=60, continue="  ")


###################################################
### code chunk number 3: ST_tutorial.Rnw:232-235
###################################################
library(SpatioTemporal)
library(plotrix) 
library(maps)


###################################################
### code chunk number 4: ST_tutorial.Rnw:247-249
###################################################
data(mesa.data.raw, package="SpatioTemporal")
str(mesa.data.raw,1)


###################################################
### code chunk number 5: ST_tutorial.Rnw:261-262
###################################################
head(mesa.data.raw$X)


###################################################
### code chunk number 6: ST_tutorial.Rnw:276-277
###################################################
mesa.data.raw$obs[1:6,1:5]


###################################################
### code chunk number 7: ST_tutorial.Rnw:290-291
###################################################
mesa.data.raw$lax.conc.1500[1:6,1:5]


###################################################
### code chunk number 8: ST_tutorial.Rnw:303-313
###################################################
##matrix of observations
obs <- mesa.data.raw$obs
##data.frame/matrix of covariates
covars <- mesa.data.raw$X
##list/3D-array with the spatio-temporal covariates
ST.list <- list(lax.conc.1500=mesa.data.raw$lax.conc.1500)

##create STdata object
mesa.data <- createSTdata(obs, covars, SpatioTemporal=ST.list, 
                          n.basis=2)


###################################################
### code chunk number 9: ST_tutorial.Rnw:326-327
###################################################
names(mesa.data)


###################################################
### code chunk number 10: ST_tutorial.Rnw:332-333
###################################################
head(mesa.data$covars)


###################################################
### code chunk number 11: figMap (eval = FALSE)
###################################################
## ###Plot the locations, see \autoref&fig:map;
## par(mfrow=c(1,1))
## plot(mesa.data$covars$long, mesa.data$covars$lat,
##      pch=24, bg=c("red","blue")[mesa.data$covars$type],
##      xlab="Longitude", ylab="Latitude")
## 
## ###Add the map of LA
## map("county", "california", col="#FFFF0055", fill=TRUE, 
##     add=TRUE)
## 
## ##Add a legend
## legend("bottomleft", c("AQS","FIXED"), pch=24, bty="n",
##        pt.bg=c("red","blue"))


###################################################
### code chunk number 12: ST_tutorial.Rnw:382-383
###################################################
###Plot the locations, see \autoref&fig:map;
par(mfrow=c(1,1))
plot(mesa.data$covars$long, mesa.data$covars$lat,
     pch=24, bg=c("red","blue")[mesa.data$covars$type],
     xlab="Longitude", ylab="Latitude")

###Add the map of LA
map("county", "california", col="#FFFF0055", fill=TRUE, 
    add=TRUE)

##Add a legend
legend("bottomleft", c("AQS","FIXED"), pch=24, bty="n",
       pt.bg=c("red","blue"))


###################################################
### code chunk number 13: ST_tutorial.Rnw:392-394
###################################################
head(mesa.data$trend)
head(mesa.data$trend.fnc)


###################################################
### code chunk number 14: ST_tutorial.Rnw:405-407
###################################################
cbind(mesa.data$trend.fnc(mesa.data$trend$date[1:5]), 
      mesa.data$trend[1:5,])


###################################################
### code chunk number 15: ST_tutorial.Rnw:418-419
###################################################
range(mesa.data$trend$date)


###################################################
### code chunk number 16: ST_tutorial.Rnw:426-427
###################################################
head(mesa.data$obs)


###################################################
### code chunk number 17: ST_tutorial.Rnw:454-456
###################################################
dim(mesa.data$SpatioTemp)
mesa.data$SpatioTemp[1:5,1:5,,drop=FALSE]


###################################################
### code chunk number 18: ST_tutorial.Rnw:480-481
###################################################
str(dimnames(mesa.data$SpatioTemp))


###################################################
### code chunk number 19: ST_tutorial.Rnw:485-487 (eval = FALSE)
###################################################
## as.character(sort(unique(c(mesa.data$obs$date,
##     mesa.data$trend$date))))


###################################################
### code chunk number 20: ST_tutorial.Rnw:491-492
###################################################
dimnames(mesa.data$SpatioTemp)[[3]]


###################################################
### code chunk number 21: ST_tutorial.Rnw:500-501
###################################################
print(mesa.data)


###################################################
### code chunk number 22: figTimeSpace (eval = FALSE)
###################################################
## ###Plot when observations occurr, see \autoref&fig:time_space;
## par(mfcol=c(1,1), mar=c(4.3,4.3,1,1))
## plot(mesa.data, "loc")


###################################################
### code chunk number 23: ST_tutorial.Rnw:532-533
###################################################
###Plot when observations occurr, see \autoref&fig:time_space;
par(mfcol=c(1,1), mar=c(4.3,4.3,1,1))
plot(mesa.data, "loc")


###################################################
### code chunk number 24: ST_tutorial.Rnw:549-551
###################################################
LUR <-  list(~log10.m.to.a1+s2000.pop.div.10000+km.to.coast,
             ~km.to.coast, ~km.to.coast)


###################################################
### code chunk number 25: ST_tutorial.Rnw:555-557
###################################################
cov.beta <- list(covf="exp", nugget=FALSE)
cov.nu <- list(covf="exp", nugget=~type)


###################################################
### code chunk number 26: ST_tutorial.Rnw:571-573
###################################################
locations <- list(coords=c("x","y"), long.lat=c("long","lat"), 
  others="type")


###################################################
### code chunk number 27: ST_tutorial.Rnw:580-586
###################################################
mesa.model <- createSTmodel(mesa.data, LUR=LUR, 
                            ST="lax.conc.1500",
                            cov.beta=cov.beta, 
                            cov.nu=cov.nu, 
                            locations=locations)
print(mesa.model)


###################################################
### code chunk number 28: ST_tutorial.Rnw:603-607
###################################################
cov.beta2 <- list(covf=c("exp","exp2","iid"), 
                  nugget=c(FALSE,FALSE,TRUE))
mesa.model2 <- updateCovf(mesa.model, cov.beta=cov.beta2)
print(mesa.model2)


###################################################
### code chunk number 29: ST_tutorial.Rnw:619-621
###################################################
model.dim <- loglikeSTdim(mesa.model)
str(model.dim)


###################################################
### code chunk number 30: ST_tutorial.Rnw:632-634
###################################################
x.init <- cbind(c( rep(2, model.dim$nparam.cov-1), 0),
                c( rep(c(1,-3), model.dim$m+1), -3, 0))


###################################################
### code chunk number 31: ST_tutorial.Rnw:650-652
###################################################
rownames(x.init) <- loglikeSTnames(mesa.model, all=FALSE)
x.init


###################################################
### code chunk number 32: ST_tutorial.Rnw:657-659 (eval = FALSE)
###################################################
## est.mesa.model <- estimate(mesa.model, x.init, 
##                            type="p", hessian.all=TRUE)


###################################################
### code chunk number 33: ST_tutorial.Rnw:662-663
###################################################
data(est.mesa.model, package="SpatioTemporal")


###################################################
### code chunk number 34: ST_tutorial.Rnw:693-695
###################################################
loglikeST(est.mesa.model$res.best$par, mesa.model)
est.mesa.model$res.best$value


###################################################
### code chunk number 35: ST_tutorial.Rnw:702-703
###################################################
print(est.mesa.model)


###################################################
### code chunk number 36: ST_tutorial.Rnw:718-719
###################################################
names(est.mesa.model)


###################################################
### code chunk number 37: ST_tutorial.Rnw:726-729
###################################################
names(est.mesa.model$res.best)
names(est.mesa.model$res.all[[1]])
names(est.mesa.model$res.all[[2]])


###################################################
### code chunk number 38: ST_tutorial.Rnw:773-774
###################################################
coef(est.mesa.model, pars="cov")[,1:2]


###################################################
### code chunk number 39: ST_tutorial.Rnw:796-798 (eval = FALSE)
###################################################
## pred.mesa.model <- predict(mesa.model, est.mesa.model, 
##                            pred.var=TRUE)


###################################################
### code chunk number 40: ST_tutorial.Rnw:801-802
###################################################
data(pred.mesa.model, package="SpatioTemporal")


###################################################
### code chunk number 41: ST_tutorial.Rnw:806-807
###################################################
names(pred.mesa.model)


###################################################
### code chunk number 42: ST_tutorial.Rnw:810-811
###################################################
print(pred.mesa.model)


###################################################
### code chunk number 43: ST_tutorial.Rnw:861-862
###################################################
beta <- estimateBetaFields(mesa.model)


###################################################
### code chunk number 44: figBetaFields (eval = FALSE)
###################################################
## par(mfrow=c(2,2), mar=c(3.3,3.3,1.5,1), mgp=c(2,1,0), pty="s")
## for(i in 1:3){
##   plotCI(x=beta$beta[,i], y=pred.mesa.model$beta$EX[,i],
##          uiw=1.96*beta$beta.sd[,i], err="x",
##          main=paste("Beta-field for f", i, "(t)", sep=""),
##          xlab="Empirical estimate", 
##          ylab="Spatio-Temporal Model",
##          pch=NA, sfrac=0.005, asp=1)
##   plotCI(x=beta$beta[,i], y=pred.mesa.model$beta$EX[,i],
##          uiw=1.96*sqrt(pred.mesa.model$beta$VX[,i]),
##          add=TRUE, pch=NA, sfrac=0.005)
##   abline(0, 1, col="grey")
## }


###################################################
### code chunk number 45: ST_tutorial.Rnw:890-891
###################################################
par(mfrow=c(2,2), mar=c(3.3,3.3,1.5,1), mgp=c(2,1,0), pty="s")
for(i in 1:3){
  plotCI(x=beta$beta[,i], y=pred.mesa.model$beta$EX[,i],
         uiw=1.96*beta$beta.sd[,i], err="x",
         main=paste("Beta-field for f", i, "(t)", sep=""),
         xlab="Empirical estimate", 
         ylab="Spatio-Temporal Model",
         pch=NA, sfrac=0.005, asp=1)
  plotCI(x=beta$beta[,i], y=pred.mesa.model$beta$EX[,i],
         uiw=1.96*sqrt(pred.mesa.model$beta$VX[,i]),
         add=TRUE, pch=NA, sfrac=0.005)
  abline(0, 1, col="grey")
}


###################################################
### code chunk number 46: figPredTS (eval = FALSE)
###################################################
## par(mfrow=c(4,1),mar=c(2.5,2.5,2,.5))
## for(i in c(1,10,17,22)){
##   plot(pred.mesa.model, ID=i, STmodel=mesa.model, 
##        col=c("black","red","grey"), lwd=1)
##   plot(pred.mesa.model, ID=i, pred.type="EX.mu", 
##        col="green", lwd=1, add=TRUE)
##   plot(pred.mesa.model, ID=i, pred.type="EX.mu.beta", 
##        col="blue", lwd=1, add=TRUE)
## }


###################################################
### code chunk number 47: ST_tutorial.Rnw:931-932
###################################################
par(mfrow=c(4,1),mar=c(2.5,2.5,2,.5))
for(i in c(1,10,17,22)){
  plot(pred.mesa.model, ID=i, STmodel=mesa.model, 
       col=c("black","red","grey"), lwd=1)
  plot(pred.mesa.model, ID=i, pred.type="EX.mu", 
       col="green", lwd=1, add=TRUE)
  plot(pred.mesa.model, ID=i, pred.type="EX.mu.beta", 
       col="blue", lwd=1, add=TRUE)
}


###################################################
### code chunk number 48: ST_tutorial.Rnw:947-949
###################################################
Ind.cv <- createCV(mesa.model, groups=10, min.dist=.1)
Ind.cv[1:10]


###################################################
### code chunk number 49: ST_tutorial.Rnw:961-962
###################################################
table(Ind.cv)


###################################################
### code chunk number 50: ST_tutorial.Rnw:975-982
###################################################
n.obs <- table(mesa.model$obs$ID)
ID1 <- names( n.obs[n.obs>270] )
ID2 <- names( n.obs[n.obs<=270] )
Ind.cv1 <- createCV(mesa.model, groups=10, 
                  subset=ID1, Icv.vector=FALSE)
Ind.cv2 <- createCV(mesa.model, groups=10, 
                  subset=ID2, Icv.vector=FALSE)


###################################################
### code chunk number 51: ST_tutorial.Rnw:985-987
###################################################
colSums(Ind.cv1)
colSums(Ind.cv2)


###################################################
### code chunk number 52: ST_tutorial.Rnw:990-992
###################################################
Ind.cv.final <- Ind.cv1 | Ind.cv2
colSums(Ind.cv.final)


###################################################
### code chunk number 53: ST_tutorial.Rnw:995-998
###################################################
table(Ind.cv)
##easier if we sort by number of observations in each group
rbind(sort(table(Ind.cv)), sort(colSums(Ind.cv.final)))


###################################################
### code chunk number 54: ST_tutorial.Rnw:1002-1004
###################################################
ID.cv <- sapply(split(mesa.model$obs$ID, Ind.cv),unique)
print(ID.cv)


###################################################
### code chunk number 55: ST_tutorial.Rnw:1009-1010
###################################################
mesa.model$D.beta[ID.cv[[10]],ID.cv[[10]]]


###################################################
### code chunk number 56: ST_tutorial.Rnw:1023-1029
###################################################
I.col <- apply(sapply(ID.cv,
                      function(x) mesa.model$locations$ID
                      %in% x), 1, 
               function(x) if(sum(x)==1) which(x) else 0)
names(I.col) <- mesa.model$locations$ID
print(I.col)


###################################################
### code chunk number 57: figMapCV (eval = FALSE)
###################################################
## par(mfrow=c(1,1))
## plot(mesa.model$locations$long,
##      mesa.model$locations$lat,
##      pch=23+floor(I.col/max(I.col)+.5), bg=I.col, 
##      xlab="Longitude", ylab="Latitude")
## map("county", "california", col="#FFFF0055", 
##     fill=TRUE, add=TRUE)


###################################################
### code chunk number 58: ST_tutorial.Rnw:1045-1046
###################################################
par(mfrow=c(1,1))
plot(mesa.model$locations$long,
     mesa.model$locations$lat,
     pch=23+floor(I.col/max(I.col)+.5), bg=I.col, 
     xlab="Longitude", ylab="Latitude")
map("county", "california", col="#FFFF0055", 
    fill=TRUE, add=TRUE)


###################################################
### code chunk number 59: ST_tutorial.Rnw:1074-1076 (eval = FALSE)
###################################################
## x.init <- coef(est.mesa.model, pars="cov")[,c("par","init")]
## est.cv.mesa <- estimateCV(mesa.model, x.init, Ind.cv)


###################################################
### code chunk number 60: ST_tutorial.Rnw:1079-1080
###################################################
data(est.cv.mesa, package="SpatioTemporal")


###################################################
### code chunk number 61: ST_tutorial.Rnw:1084-1085
###################################################
print(est.cv.mesa)


###################################################
### code chunk number 62: ST_tutorial.Rnw:1095-1096
###################################################
head( coef(est.cv.mesa) )


###################################################
### code chunk number 63: ST_tutorial.Rnw:1117-1118 (eval = FALSE)
###################################################
## pred.cv.mesa <- predictCV(mesa.model, est.cv.mesa, LTA=TRUE)


###################################################
### code chunk number 64: ST_tutorial.Rnw:1121-1122
###################################################
data(pred.cv.mesa, package="SpatioTemporal")


###################################################
### code chunk number 65: ST_tutorial.Rnw:1127-1129
###################################################
print(pred.cv.mesa)
names(pred.cv.mesa)


###################################################
### code chunk number 66: ST_tutorial.Rnw:1133-1134
###################################################
str(pred.cv.mesa$pred.obs)


###################################################
### code chunk number 67: ST_tutorial.Rnw:1138-1139
###################################################
str(pred.cv.mesa$pred.all,1)


###################################################
### code chunk number 68: ST_tutorial.Rnw:1143-1144
###################################################
names(pred.mesa.model)


###################################################
### code chunk number 69: ST_tutorial.Rnw:1148-1149
###################################################
str(pred.cv.mesa$pred.LTA)


###################################################
### code chunk number 70: ST_tutorial.Rnw:1153-1154
###################################################
summary(pred.cv.mesa)


###################################################
### code chunk number 71: ST_tutorial.Rnw:1172-1175
###################################################
I.season <- as.factor(as.POSIXlt(pred.cv.mesa$pred.obs$date)$mon+1)
levels(I.season) <- c(rep("Winter",2), rep("Spring",3), 
                      rep("Summer",3), rep("Fall",3), "Winter") 


###################################################
### code chunk number 72: ST_tutorial.Rnw:1179-1180
###################################################
table(I.season)


###################################################
### code chunk number 73: figPredCvQQplot (eval = FALSE)
###################################################
## par(mfrow=c(1,2), mar=c(3,2,1,1), pty="s")
## qqnorm(pred.cv.mesa, col=I.season, line=2)
## qqnorm(pred.cv.mesa, norm=TRUE, main="Normalised residuals",
##        col=I.season)
## legend("bottomright", legend=as.character(levels(I.season)),
##        pch=1, col=1:nlevels(I.season))


###################################################
### code chunk number 74: ST_tutorial.Rnw:1202-1203
###################################################
par(mfrow=c(1,2), mar=c(3,2,1,1), pty="s")
qqnorm(pred.cv.mesa, col=I.season, line=2)
qqnorm(pred.cv.mesa, norm=TRUE, main="Normalised residuals",
       col=I.season)
legend("bottomright", legend=as.character(levels(I.season)),
       pch=1, col=1:nlevels(I.season))


###################################################
### code chunk number 75: figCVresids (eval = FALSE)
###################################################
## par(mfcol=c(2,1),mar=c(4.5,4.5,2,2))
## scatterPlot(pred.cv.mesa, trend=1, group=I.season, col=c(2:5,1), 
##             xlab="First temporal smooth", type="res",
##             STdata=mesa.model, main="CV Residuals - All data")
## scatterPlot(pred.cv.mesa, covar="log10.m.to.a1", group=I.season, 
##             col=c(2:5,1), STdata=mesa.model, type="res",
##             main="CV Residuals - All data")
## legend("topleft", levels(I.season), col=c(2:5), pch=1, cex=.75)


###################################################
### code chunk number 76: figCVresidsByType (eval = FALSE)
###################################################
## par(mfcol=c(1,2),mar=c(4.5,4.5,2,2))
## scatterPlot(pred.cv.mesa, covar="log10.m.to.a1", group=I.season, 
##             subset=with(mesa.data$covars, ID[type=="AQS"]),
##             col=c(2:5,1), lty=c(rep(2,4),1), type="res",
##             STdata=mesa.model, main="AQS sites")
## legend("topleft", levels(I.season), col=c(2:5), pch=1, cex=.75)
## ##and for the FIXED sites
## scatterPlot(pred.cv.mesa, covar="log10.m.to.a1", group=I.season, 
##             subset=with(mesa.data$covars, ID[type=="FIXED"]),
##             col=c(2:5,1), lty=c(rep(2,4),1), type="res",
##             STdata=mesa.model, main="FIXED sites")


###################################################
### code chunk number 77: ST_tutorial.Rnw:1255-1256
###################################################
par(mfcol=c(2,1),mar=c(4.5,4.5,2,2))
scatterPlot(pred.cv.mesa, trend=1, group=I.season, col=c(2:5,1), 
            xlab="First temporal smooth", type="res",
            STdata=mesa.model, main="CV Residuals - All data")
scatterPlot(pred.cv.mesa, covar="log10.m.to.a1", group=I.season, 
            col=c(2:5,1), STdata=mesa.model, type="res",
            main="CV Residuals - All data")
legend("topleft", levels(I.season), col=c(2:5), pch=1, cex=.75)


###################################################
### code chunk number 78: ST_tutorial.Rnw:1263-1264
###################################################
par(mfcol=c(1,2),mar=c(4.5,4.5,2,2))
scatterPlot(pred.cv.mesa, covar="log10.m.to.a1", group=I.season, 
            subset=with(mesa.data$covars, ID[type=="AQS"]),
            col=c(2:5,1), lty=c(rep(2,4),1), type="res",
            STdata=mesa.model, main="AQS sites")
legend("topleft", levels(I.season), col=c(2:5), pch=1, cex=.75)
##and for the FIXED sites
scatterPlot(pred.cv.mesa, covar="log10.m.to.a1", group=I.season, 
            subset=with(mesa.data$covars, ID[type=="FIXED"]),
            col=c(2:5,1), lty=c(rep(2,4),1), type="res",
            STdata=mesa.model, main="FIXED sites")


###################################################
### code chunk number 79: ST_tutorial.Rnw:1302-1304
###################################################
##clean up
rm(list=ls())


###################################################
### code chunk number 80: ST_tutorial.Rnw:1312-1320
###################################################
##libraries
library(SpatioTemporal)
library(plotrix) 

##load data
data(mesa.data.raw)
data(mesa.model)
data(est.mesa.model)


###################################################
### code chunk number 81: ST_tutorial.Rnw:1326-1331
###################################################
mesa.data.raw$obs <- 
  mesa.data.raw$obs[,!(colnames(mesa.data.raw$obs) %in%
                       c("60595001", "LC003"))] 
mesa.data <- with(mesa.data.raw, 
                  createSTdata(obs, X, n.basis=2))


###################################################
### code chunk number 82: ST_tutorial.Rnw:1335-1338
###################################################
mesa.data.org <- mesa.data
T <- with(mesa.data$trend, seq(min(date), max(date), by=7))
mesa.data <- updateTrend(mesa.data, n.basis=2, extra.dates=T)


###################################################
### code chunk number 83: ST_tutorial.Rnw:1344-1346
###################################################
print(mesa.data)
print(mesa.data.org)


###################################################
### code chunk number 84: ST_tutorial.Rnw:1349-1358
###################################################
par(mfrow=c(2,1), mar=c(2,2,2,1))
plot(mesa.data$trend$date, mesa.data$trend$V1, 
     xlab="", ylab="", main="Trend 1")
points(mesa.data.org$trend$date, mesa.data.org$trend$V1, 
       col="red", pch=3)
plot(mesa.data$trend$date, mesa.data$trend$V2,
     xlab="", ylab="", main="Trend 2")
points(mesa.data.org$trend$date, mesa.data.org$trend$V2,
       col="red", pch=3)


###################################################
### code chunk number 85: ST_tutorial.Rnw:1365-1373
###################################################
mesa.model.1 <- createSTmodel(mesa.data,
  LUR=mesa.model$LUR.list, cov.beta=mesa.model$cov.beta, 
  cov.nu=mesa.model$cov.nu, 
  locations=mesa.model$locations.list)
mesa.model.2 <- createSTmodel(mesa.data.org, 
  LUR=mesa.model$LUR.list, cov.beta=mesa.model$cov.beta, 
  cov.nu=mesa.model$cov.nu, 
  locations=mesa.model$locations.list, strip=TRUE)


###################################################
### code chunk number 86: ST_tutorial.Rnw:1376-1378
###################################################
print(mesa.model.1)
print(mesa.model.2)


###################################################
### code chunk number 87: ST_tutorial.Rnw:1381-1383
###################################################
str(loglikeSTdim(mesa.model.1))
str(loglikeSTdim(mesa.model.2))


###################################################
### code chunk number 88: ST_tutorial.Rnw:1386-1388
###################################################
dim(mesa.model.1$trend)
dim(mesa.model.2$trend)


###################################################
### code chunk number 89: ST_tutorial.Rnw:1399-1400
###################################################
x <- coef(est.mesa.model, pars="cov")$par


###################################################
### code chunk number 90: ST_tutorial.Rnw:1404-1406
###################################################
E.1 <- predict(mesa.model.1, x, pred.var=FALSE)
E.2 <- predict(mesa.model.2, x, pred.var=FALSE)


###################################################
### code chunk number 91: ST_tutorial.Rnw:1411-1415
###################################################
colnames(E.1$EX)
str(E.1$EX)
colnames(E.2$EX)
str(E.2$EX)


###################################################
### code chunk number 92: ST_tutorial.Rnw:1418-1419
###################################################
range(E.1$EX[rownames(E.2$EX),colnames(E.2$EX)] - E.2$EX)


###################################################
### code chunk number 93: ST_tutorial.Rnw:1423-1424
###################################################
E.3 <- predict(mesa.model.2, x, STdata=mesa.data, pred.var=FALSE)


###################################################
### code chunk number 94: ST_tutorial.Rnw:1430-1433
###################################################
colnames(E.3$EX)
str(E.3$EX)
all.equal(E.3,E.1)


###################################################
### code chunk number 95: ST_tutorial.Rnw:1435-1436
###################################################
stopifnot( all.equal(E.3,E.1) )


###################################################
### code chunk number 96: ST_tutorial.Rnw:1446-1449
###################################################
LTA <- with(mesa.model.1$trend, split(date,as.POSIXlt(date)$year+1900))
str(LTA)
lapply(LTA[1:3], range)


###################################################
### code chunk number 97: ST_tutorial.Rnw:1453-1456
###################################################
ID <- mesa.model.1$locations$ID
LTA <- rep(list(LTA), length(ID))
names(LTA) <- ID


###################################################
### code chunk number 98: ST_tutorial.Rnw:1460-1461
###################################################
E.1.LTA <- predict(mesa.model.1, x, pred.var=FALSE, LTA=LTA)


###################################################
### code chunk number 99: ST_tutorial.Rnw:1465-1466
###################################################
head(E.1.LTA$LTA)


###################################################
### code chunk number 100: ST_tutorial.Rnw:1477-1479
###################################################
E.1.LTA.alt <- sapply(split(E.1$EX[,1], as.POSIXlt(rownames(E.1$EX))$year),mean)
cbind(E.1.LTA.alt, E.1.LTA$LTA$EX[1:11])


###################################################
### code chunk number 101: ST_tutorial.Rnw:1481-1482
###################################################
stopifnot( max(abs(E.1.LTA.alt- E.1.LTA$LTA$EX[1:11]))<1e-14 )


###################################################
### code chunk number 102: ST_tutorial.Rnw:1486-1488
###################################################
##clean up
rm(list=ls())


###################################################
### code chunk number 103: ST_tutorial.Rnw:1499-1507
###################################################
##Load libraries
library(SpatioTemporal)
library(plotrix) 
library(maps)

##and data
data(mesa.model)
data(est.mesa.model)


###################################################
### code chunk number 104: ST_tutorial.Rnw:1513-1515
###################################################
x <- coef(est.mesa.model)$par
sim.data <- simulate(mesa.model, nsim=4, x=x)


###################################################
### code chunk number 105: ST_tutorial.Rnw:1518-1520
###################################################
names(sim.data)
str(sim.data,1)


###################################################
### code chunk number 106: ST_tutorial.Rnw:1527-1534
###################################################
mesa.data.sim <- list()
for(i in 1:length(sim.data$obs)){
  ##copy the mesa.data.model object
  mesa.data.sim[[i]] <- mesa.model
  ##replace observations with the simulated data
  mesa.data.sim[[i]]$obs <- sim.data$obs[[i]]
}


###################################################
### code chunk number 107: ST_tutorial.Rnw:1540-1545
###################################################
data(mesa.data.raw)
mesa.data.raw$X <- mesa.data.raw$X[mesa.data.raw$X[,"ID"]=="60590001",]
mesa.data <- createSTdata(obs=NULL, covars=mesa.data.raw$X,
   extra.dates=as.Date(mesa.model$trend$date),
   SpatioTemporal=list(lax.conc.1500=mesa.data.raw$lax.conc.1500))


###################################################
### code chunk number 108: ST_tutorial.Rnw:1551-1555
###################################################
E <- list()
for(i in 1:length(sim.data$obs)){
  E[[i]] <- predict(mesa.data.sim[[i]], x, STdata=mesa.data)
}


###################################################
### code chunk number 109: ST_tutorial.Rnw:1564-1573
###################################################
par(mfrow=c(2,2),mar=c(2.5,2.5,2,.5))
for(i in 1:4){
  ##plot predictions, but not the observations
  plot(E[[i]])
  ##add the simulated data (i.e. observations + 
  ##simulated values at points where we've predicted)
  lines(as.Date(rownames(sim.data$X)), 
        sim.data$X[,mesa.data$covars$ID,i], col="red")
}


###################################################
### code chunk number 110: ST_tutorial.Rnw:1577-1579
###################################################
##clean up
rm(list=ls())


###################################################
### code chunk number 111: ST_tutorial.Rnw:1588-1594
###################################################
library(SpatioTemporal)
library(plotrix) 

##load data
data(mesa.model)
data(est.mesa.model)


###################################################
### code chunk number 112: ST_tutorial.Rnw:1605-1609
###################################################
##parameters
x <- coef(est.mesa.model)
##and Hessian
H <- est.mesa.model$res.best$hessian.all


###################################################
### code chunk number 113: ST_tutorial.Rnw:1612-1614 (eval = FALSE)
###################################################
## MCMC.mesa.model <- MCMC(mesa.model, x$par, N = 2500, 
##                         Hessian.prop = H)


###################################################
### code chunk number 114: ST_tutorial.Rnw:1617-1618
###################################################
data(MCMC.mesa.model)


###################################################
### code chunk number 115: ST_tutorial.Rnw:1624-1626
###################################################
print(MCMC.mesa.model)
names(MCMC.mesa.model)


###################################################
### code chunk number 116: ST_tutorial.Rnw:1629-1630
###################################################
summary(MCMC.mesa.model)


###################################################
### code chunk number 117: ST_tutorial.Rnw:1636-1640
###################################################
par(mfrow=c(4,1),mar=c(2,2,2.5,.5))
for(i in c(4,9,13,15)){
  plot(MCMC.mesa.model, i, ylab="", xlab="", type="l")
}


###################################################
### code chunk number 118: ST_tutorial.Rnw:1646-1653
###################################################
dens <- density(MCMC.mesa.model, estSTmodel=x)

##plots for all covariance parameters
par(mfrow=c(3,3),mar=c(4,4,2.5,.5))
for(i in 9:17){
  plot(dens, i, norm.col="red")
}


