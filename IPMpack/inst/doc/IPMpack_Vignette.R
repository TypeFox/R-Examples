### R code from vignette source 'IPMpack_Vignette.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: IPMpack_Vignette.Rnw:157-159
###################################################
library(IPMpack)



###################################################
### code chunk number 2: IPMpack_Vignette.Rnw:161-162
###################################################
options(continue="  ") # removes those stupid plus signs from R continuations


###################################################
### code chunk number 3: IPMpack_Vignette.Rnw:164-171
###################################################
library(IPMpack)
#source("/Applications/eclipse//IPMpack/R/IPMpack-Util.r")
#source("/Applications/eclipse//IPMpack/R/IPMpack-Analyses.r")
#source("/Applications/eclipse//IPMpack/R/IPMpack-Vital_Rate_Obj.r")
#source("/Applications/eclipse//IPMpack/R/IPMpack-Hidden.r")
#source("/Applications/eclipse//IPMpack/R/IPMpack-Classes_and_Methods.r")
#source("/Applications/eclipse//IPMpack/R/IPMpack-Matrices.r")


###################################################
### code chunk number 4: IPMpack_Vignette.Rnw:175-176
###################################################
dff <- generateData()


###################################################
### code chunk number 5: IPMpack_Vignette.Rnw:179-180
###################################################
head(dff)


###################################################
### code chunk number 6: fig0
###################################################
plot(dff$size, dff$sizeNext, xlab = "Size at t", ylab = "Size at t+1")


###################################################
### code chunk number 7: fig0
###################################################
plot(dff$size, dff$sizeNext, xlab = "Size at t", ylab = "Size at t+1")


###################################################
### code chunk number 8: IPMpack_Vignette.Rnw:215-217
###################################################
gr1 <- makeGrowthObj(dataf = dff, Formula = sizeNext~size+size2)
sv1 <- makeSurvObj(dff, Formula = surv~size+size2)


###################################################
### code chunk number 9: IPMpack_Vignette.Rnw:222-223
###################################################
gr1


###################################################
### code chunk number 10: fig1
###################################################
par(mfrow = c(1, 2), bty = "l", pty = "m")
p1 <- picGrow(dff, gr1)
p2 <- picSurv(dff, sv1, ncuts = 30)


###################################################
### code chunk number 11: fig1
###################################################
par(mfrow = c(1, 2), bty = "l", pty = "m")
p1 <- picGrow(dff, gr1)
p2 <- picSurv(dff, sv1, ncuts = 30)


###################################################
### code chunk number 12: IPMpack_Vignette.Rnw:241-245
###################################################
Pmatrix <- makeIPMPmatrix(nBigMatrix = 50, 
                            minSize = -5, maxSize = 35, 
                            growObj = gr1, survObj = sv1, 
                            correction = "constant")


###################################################
### code chunk number 13: triDiag
###################################################
diagnosticsPmatrix(Pmatrix, gr1, sv1, dff,correction = "constant")


###################################################
### code chunk number 14: IPMpack_Vignette.Rnw:285-286
###################################################
slotNames(Pmatrix)


###################################################
### code chunk number 15: IPMpack_Vignette.Rnw:289-290
###################################################
Pmatrix@meshpoints


###################################################
### code chunk number 16: fig2
###################################################
par(mfrow=c(1,2),bty="l",pty="s")
persp(Pmatrix@meshpoints, Pmatrix@meshpoints, Pmatrix, ticktype = "detailed", 
  	  xlab = "Size at t", ylab = "Size at t+1", zlab="")
image(Pmatrix@meshpoints, Pmatrix@meshpoints, t(Pmatrix), 
		  xlab = "Size at t", ylab = "Size at t+1")
contour(Pmatrix@meshpoints, Pmatrix@meshpoints, t(Pmatrix), add=TRUE)


###################################################
### code chunk number 17: IPMpack_Vignette.Rnw:310-312
###################################################
LE <- meanLifeExpect(Pmatrix)
pTime <- passageTime(mean(dff$size, na.rm = TRUE), Pmatrix)


###################################################
### code chunk number 18: fig3
###################################################
par(mfrow = c(1, 2), bty = "l")
plot(Pmatrix@meshpoints, LE, type = "l", xlab = "Size", 
     ylab = "Mean life expectancy", 
     xlim = range(dff$size, na.rm = TRUE), 
     ylim = range(LE[Pmatrix@meshpoints<max(dff$size, na.rm = TRUE)]))   
plot(Pmatrix@meshpoints, pTime, type = "l", 
     xlab = "Size at start", ylab = "Time to reach chosen size", 
     xlim = range(Pmatrix@meshpoints[Pmatrix@meshpoints<mean(dff$size, na.rm = TRUE)]), 
     ylim = range(pTime[Pmatrix@meshpoints<max(dff$size, na.rm = TRUE)]))
abline(v = mean(dff$size, na.rm = TRUE), col = 2)  #show the target size in red


###################################################
### code chunk number 19: IPMpack_Vignette.Rnw:337-344
###################################################
fv1 <- makeFecObj(dff, Formula = fec~size, 	
                  Family = "gaussian", 
                  Transform = "log")
Fmatrix <- makeIPMFmatrix(nBigMatrix = 50, minSize = -5,
                            maxSize = 35, 
                            fecObj = fv1, 
                            correction = "constant")


###################################################
### code chunk number 20: IPMpack_Vignette.Rnw:360-364
###################################################
IPM <- Pmatrix + Fmatrix
Re(eigen(IPM)$value[1])
sensitivity <- sens(IPM)
elasticity <- elas(IPM)


###################################################
### code chunk number 21: fig4
###################################################
par(mfrow = c(1, 3), bty = "l", pty = "s")
plot(Pmatrix@meshpoints, abs(Re(eigen(IPM)$vector[ , 1])), 
		 type = "l", xlab = "Size", 
     ylab = "Stable size structure")
image(Pmatrix@meshpoints, Pmatrix@meshpoints, log(sensitivity), 
		 xlab = "Size at t", 
      ylab = "Size at t+1", main = "Sensitivity")
image(Pmatrix@meshpoints, Pmatrix@meshpoints, log(elasticity), 
		 xlab = "Size at t", 
      ylab = "Size at t+1", main = "Elasticity")


###################################################
### code chunk number 22: IPMpack_Vignette.Rnw:386-389
###################################################
res <- sensParams(growObj = gr1, survObj = sv1, fecObj = fv1, 
                  nBigMatrix = 50, minSize = -5, maxSize = 15)
res


###################################################
### code chunk number 23: fig4a
###################################################
par(mfrow = c(2, 1), bty = "l", pty = "m")
barplot(res$sens, main = expression("Parameter sensitivity of "*lambda), 
		    las = 2, cex.names = 0.5)
barplot(res$elas, main = expression("Parameter elasticity of "*lambda), 
		    las = 2, cex.names = 0.5)


###################################################
### code chunk number 24: fig4a
###################################################
par(mfrow = c(2, 1), bty = "l", pty = "m")
barplot(res$sens, main = expression("Parameter sensitivity of "*lambda), 
		    las = 2, cex.names = 0.5)
barplot(res$elas, main = expression("Parameter elasticity of "*lambda), 
		    las = 2, cex.names = 0.5)


###################################################
### code chunk number 25: IPMpack_Vignette.Rnw:418-419
###################################################
dff <- generateData(type="discrete")


###################################################
### code chunk number 26: IPMpack_Vignette.Rnw:422-423
###################################################
table(dff$stage)


###################################################
### code chunk number 27: IPMpack_Vignette.Rnw:431-436
###################################################
fv1 <- makeFecObj(dataf = dff, Transform = "log", 
                  offspringSplitter = data.frame(continuous = 0.2, 
                  dormant = 0, seedAge1 = 0.8, seedOld = 0), 
                  fecByDiscrete = data.frame(dormant = 0, 
                  seedAge1 = 0, seedOld = 0))


###################################################
### code chunk number 28: IPMpack_Vignette.Rnw:445-449
###################################################
Fmatrix <- makeIPMFmatrix(fecObj = fv1, nBigMatrix = 5, 
                            minSize = min(dff$size, na.rm = TRUE), 
                            maxSize = max(dff$size, na.rm = TRUE), 
                            correction = "constant")


###################################################
### code chunk number 29: IPMpack_Vignette.Rnw:453-456
###################################################
gr1 <- makeGrowthObj(dataf = dff, 
                      Formula = sizeNext~size)
sv1 <- makeSurvObj(dff,  Formula = surv~size)


###################################################
### code chunk number 30: IPMpack_Vignette.Rnw:459-460
###################################################
discTrans <- makeDiscreteTrans(dff)


###################################################
### code chunk number 31: IPMpack_Vignette.Rnw:466-473
###################################################
Pmatrix <- makeIPMPmatrix(nBigMatrix = 5, 	
                            minSize = min(dff$size, na.rm = TRUE), 
                            maxSize = max(dff$size, na.rm = TRUE), 
                            growObj = makeGrowthObj(dff), 
                            survObj = makeSurvObj(dff), 
                            discreteTrans = discTrans, 
                            correction = "constant")


###################################################
### code chunk number 32: IPMpack_Vignette.Rnw:476-478
###################################################
print(Pmatrix)
print(Fmatrix)


###################################################
### code chunk number 33: IPMpack_Vignette.Rnw:482-483
###################################################
print(Pmatrix+Fmatrix)


###################################################
### code chunk number 34: IPMpack_Vignette.Rnw:497-500
###################################################
dff <- generateData()
env1 <- makeEnvObj(dff) 
env1


###################################################
### code chunk number 35: IPMpack_Vignette.Rnw:504-506
###################################################
gr1 <- makeGrowthObj(dff, Formula = sizeNext~size+covariate)
sv1 <- makeSurvObj(dff, Formula = surv~size+covariate)


###################################################
### code chunk number 36: IPMpack_Vignette.Rnw:518-523
###################################################
Pmatrix <- makeCompoundPmatrix(nBigMatrix = 50, minSize = -5, 
                                 maxSize = 35, 
                                 envMatrix = env1, growObj = gr1, 
                                 survObj = sv1, 
                                 correction = "constant")


###################################################
### code chunk number 37: figCompound
###################################################
image(1:nrow(Pmatrix), 1:ncol(Pmatrix), t(log(Pmatrix)), 
	xlab = "Continuous stage (e.g. size) at t", 
		ylab = "Continuous stage (e.g. size) at t+1", axes = FALSE)
axis(1, at = 1:nrow(Pmatrix), lab = round(rep(Pmatrix@meshpoints,Pmatrix@nEnvClass), 2))
axis(2, at = 1:nrow(Pmatrix), lab = round(rep(Pmatrix@meshpoints,Pmatrix@nEnvClass), 2))
abline(h = length(Pmatrix@meshpoints) * (1:Pmatrix@nEnvClass))
abline(v = length(Pmatrix@meshpoints) * (1:Pmatrix@nEnvClass))


###################################################
### code chunk number 38: IPMpack_Vignette.Rnw:547-548
###################################################
pTimes <- stochPassageTime(Pmatrix@meshpoints[15], Pmatrix, env1)


###################################################
### code chunk number 39: fig5
###################################################
par(mfrow = c(1, 1), bty = "l")
#plot(Pmatrix@meshpoints, LEs[1, ], type = "l", xlab = "Size", 
#     ylab = "Mean life expectancy", ylim = c(0, 8), 
#     xlim = range(dff$size, na.rm = TRUE))   
#for (k in 1:Pmatrix@nEnvClass)  {
#	points(Pmatrix@meshpoints, LEs[k, ], type = "l", col = k)
#	}
plot(Pmatrix@meshpoints, pTimes[1:Pmatrix@nBigMatrix], 
     type  = "l", xlab = "Current size", 
     ylab = "Time to reach chosen size", 
     xlim = range(Pmatrix@meshpoints[1:14]), ylim = c(0, 8))
 for (i in 1:Pmatrix@nEnvClass) {
	points(Pmatrix@meshpoints, 
               pTimes[((i-1) * Pmatrix@nBigMatrix+1):(Pmatrix@nBigMatrix*i)], 
               type = "l", col = i)
  }


###################################################
### code chunk number 40: sampleSequentialIPMs
###################################################
IPMlist <- sampleSequentialIPMs(dataf = dff, nBigMatrix = 25, minSize = -5, 
                        maxSize = 35, explSurv = surv~size+covariate, 
                        explGrow = sizeNext~size+size2+covariate, 
                        explFec = fec~size, Transform="log",correction = "constant")


###################################################
### code chunk number 41: stochGrowth
###################################################
stochGrowthRateSampleList(listIPMmatrix = IPMlist, 
                          nRunIn = 30, tMax = 50)


###################################################
### code chunk number 42: IPMpack_Vignette.Rnw:611-613
###################################################
dff <- generateData()
dff$covariate <- sample(c(1:4),nrow(dff),replace=TRUE)


###################################################
### code chunk number 43: IPMpack_Vignette.Rnw:619-624
###################################################
gr1 <- makeGrowthObj(dataf = dff,
                     Formula=sizeNext~size+covariate,
                     regType="constantVar",
                     Family="gaussian")
sv1 <- makeSurvObj(dff, Formula = surv~size+covariate)


###################################################
### code chunk number 44: IPMpack_Vignette.Rnw:628-633
###################################################
n.age.classes <- max(dff$covariate,na.rm=TRUE)
ageMat <- new("envMatrix", nEnvClass = n.age.classes)
ageMat@.Data <- matrix(0,n.age.classes,n.age.classes)
ageMat@.Data[cbind(2:n.age.classes,1:(n.age.classes-1))] <- 1
ageMat@.Data[n.age.classes,n.age.classes] <- 1


###################################################
### code chunk number 45: IPMpack_Vignette.Rnw:638-639
###################################################
print(ageMat)


###################################################
### code chunk number 46: IPMpack_Vignette.Rnw:644-651
###################################################
Pmatrix <- makeCompoundPmatrix(nBigMatrix = 50, minSize = -5,
                               maxSize = 35,
                               envMatrix = ageMat,
							   nEnvClass = 4,
                               growObj = gr1,
                               survObj = sv1,
                               correction = "constant")


###################################################
### code chunk number 47: figCompound
###################################################
image(Pmatrix[,],
	xlab = "Continuous stage (e.g. size) at t",
		ylab = "Continuous stage (e.g. size) at t+1", axes = FALSE)


###################################################
### code chunk number 48: IPMpack_Vignette.Rnw:674-677
###################################################
fv1 <- makeFecObj(dff, Formula = fec~size,
                    Family = "gaussian",
                    Transform = "log")


###################################################
### code chunk number 49: IPMpack_Vignette.Rnw:680-684
###################################################
n.age.classes <- max(dff$covariate,na.rm=TRUE)
ageMat1 <- new("envMatrix", nEnvClass = n.age.classes)
ageMat1@.Data <- matrix(0,n.age.classes,n.age.classes)
ageMat1@.Data[1,2:n.age.classes] <- 1


###################################################
### code chunk number 50: IPMpack_Vignette.Rnw:687-693
###################################################
Fmatrix <- makeCompoundFmatrix(nBigMatrix = 50, minSize = -5,
                               maxSize = 35,
                               envMatrix = ageMat1,
							   nEnvClass = 4,
                               fecObj = fv1,
                               correction = "constant")


###################################################
### code chunk number 51: figCompound
###################################################
IPM <- Pmatrix+Fmatrix
image(log(IPM),
	xlab = "Continuous stage (e.g. size) at t",
		ylab = "Continuous stage (e.g. size) at t+1", axes = FALSE)


###################################################
### code chunk number 52: buildNew
###################################################
dff <- generateData(type="stochastic")
sv1 <- makeSurvObj(dataf = dff, 
                          Formula = surv~size+covariate1+covariate3)
gr1 <- makeGrowthObj(dataf = dff, 
                             Formula = sizeNext~size+covariate1+covariate2)
fv1 <- makeFecObj(dataf = dff, fecConstants = data.frame(1.8), 
                  Formula = fec~size, Transform = "log")


###################################################
### code chunk number 53: headDff
###################################################
head(dff)


###################################################
### code chunk number 54: grobj
###################################################
gr1


###################################################
### code chunk number 55: covTests
###################################################
tVals <- seq(1, 4, by = 1/12)
covTest <- c(1 + 0.5*sin(2*pi*tVals))
covMatTest <- data.frame(covariate1 = rnorm(length(covTest), covTest, 0.5) - 1, 
                         covariate2 = rnorm(length(covTest), covTest, 0.5) - 1, 
                         covariate3 = rnorm(length(covTest), covTest, 0.5) - 1)


###################################################
### code chunk number 56: stochGROWTH
###################################################
r <- stochGrowthRateManyCov(covariate = covMatTest, nRunIn = 12*1, 
                            tMax = length(tVals), growthObj = gr1, 
                            survObj = sv1, fecObj = fv1, nBigMatrix = 20, 
                            minSize = 2*min(dff$size, na.rm = TRUE), 
                            maxSize = 1.6*max(dff$size, na.rm = TRUE), 
                            nMicrosites = 50, correction = "constant")
print(r)


###################################################
### code chunk number 57: setClass
###################################################
setClass("growthObjSaturate", representation(paras = "numeric", sd = "numeric"))


###################################################
### code chunk number 58: fsat
###################################################
fSaturate <- function(size, pars) { 
    u <- exp(pmin(pars[1] + pars[2] * size, 50))
    u <- pars[3] * 1/(1+u)
    return(u)
}


###################################################
### code chunk number 59: wrapSat
###################################################
wrapSaturate <- function(par, dataf) { 
    pred <- fSaturate(dataf$size, par[1:3])
    ss <- sum((pred - dataf$sizeNext)^2, na.rm = TRUE)
    return(ss)
    }
tmp <- optim(c(1, 1, 1), wrapSaturate, dataf = dff, method = "Nelder-Mead")
tmp    


###################################################
### code chunk number 60: residsSd
###################################################
resids <- fSaturate(dff$size, tmp$par) - dff$sizeNext
sdSaturate <- sd(resids, na.rm = TRUE)


###################################################
### code chunk number 61: gr1new
###################################################
gr1 <- new("growthObjSaturate")
gr1@paras <- tmp$par
gr1@sd <- sdSaturate


###################################################
### code chunk number 62: newmethod
###################################################
setMethod("growth", c("numeric", "numeric", "numeric", "growthObjSaturate"), 
          function(size, sizeNext, cov, growthObj){
              mux <- fSaturate(size, growthObj@paras)
              sigmax <- growthObj@sd
              u <- dnorm(sizeNext, mux, sigmax, log = F)  
              return(u);
          })


