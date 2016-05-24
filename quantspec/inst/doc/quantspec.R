### R code from vignette source 'quantspec.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: quantspec.Rnw:19-20
###################################################
set.seed(2581)


###################################################
### code chunk number 2: quantspec.Rnw:735-736
###################################################
library("quantspec")


###################################################
### code chunk number 3: quantspec.Rnw:742-743 (eval = FALSE)
###################################################
## help("quantspec")


###################################################
### code chunk number 4: quantspec.Rnw:757-759
###################################################
Y  <- rnorm(8)
bn <- qRegEstimator(Y, levels = c(0.25, 0.5, 0.75))


###################################################
### code chunk number 5: quantspec.Rnw:766-767
###################################################
bn


###################################################
### code chunk number 6: quantspec.Rnw:779-781
###################################################
getFrequencies(bn)
getParallel(bn)


###################################################
### code chunk number 7: quantspec.Rnw:790-791
###################################################
getValues(bn, levels = c(0.25, 0.5))


###################################################
### code chunk number 8: quantspec.Rnw:799-800 (eval = FALSE)
###################################################
## help("getValues-FreqRep")


###################################################
### code chunk number 9: plotbn (eval = FALSE)
###################################################
## dn <- clippedFT(rnorm(32), levels = seq(0.05, 0.95, 0.05))
## plot(dn, frequencies = 2 * pi * (0:64) / 32, levels = c(0.25, 0.5))


###################################################
### code chunk number 10: quantspec.Rnw:820-821
###################################################
dn <- clippedFT(rnorm(32), levels = seq(0.05, 0.95, 0.05))
plot(dn, frequencies = 2 * pi * (0:64) / 32, levels = c(0.25, 0.5))


###################################################
### code chunk number 11: quantspec.Rnw:857-860 (eval = FALSE)
###################################################
## demo("sp500")
## demo("wheatprices")
## demo("qar-simulation")


###################################################
### code chunk number 12: plotSP (eval = FALSE)
###################################################
## library("zoo")
## plot(sp500,            xlab = "time t", ylab = "", main = "")
## acf(coredata(sp500),   xlab = "lag k",  ylab = "", main = "")
## acf(coredata(sp500)^2, xlab = "lag k",  ylab = "", main = "")


###################################################
### code chunk number 13: quantspec.Rnw:909-914
###################################################
def.par <- par(no.readonly = TRUE) # save default, for resetting...
layout(matrix(c(1,1,2,3), nrow=1))
par(mar=c(5,2,2,1))
library("zoo")
plot(sp500,            xlab = "time t", ylab = "", main = "")
acf(coredata(sp500),   xlab = "lag k",  ylab = "", main = "")
acf(coredata(sp500)^2, xlab = "lag k",  ylab = "", main = "")
par(def.par)  #- reset to default


###################################################
### code chunk number 14: SPquantilePG (eval = FALSE)
###################################################
## CR <- quantilePG(sp500, levels.1 = c(0.05, 0.5, 0.95),
##   type = "clipped", type.boot = "mbb", B = 250, l = 32)
## freq <- getFrequencies(CR)
## plot(CR, levels = c(0.05, 0.5, 0.95),
##   frequencies = freq[freq > 0 & freq <= pi],
##   ylab = expression({I[list(n, R)]^{list(tau[1], tau[2])}}(omega)))


###################################################
### code chunk number 15: SPquantilePG2 (eval = FALSE)
###################################################
## plot(CR, levels = c(0.05, 0.5, 0.95),
##   frequencies = freq[freq > 0 & freq <= pi/5],
##   ylab = expression({I[list(n, R)]^{list(tau[1], tau[2])}}(omega)))


###################################################
### code chunk number 16: SPsmoothedPG (eval = FALSE)
###################################################
## sPG <- smoothedPG(CR, weight = kernelWeight(W = W1, bw = 0.07))
## plot(sPG, levels = c(0.05, 0.5, 0.95), type.scaling = "individual",
##   frequencies = freq[freq > 0 & freq <= pi], ptw.CIs = 0.1,
##   ylab = expression(hat(G)[list(n, R)](list(tau[1], tau[2], omega))))


###################################################
### code chunk number 17: SPsmoothedPGboot (eval = FALSE)
###################################################
## plot(sPG, levels = c(0.05, 0.5, 0.95), type.scaling = "real-imaginary",
##   ptw.CIs = 0.1, type.CIs = "boot.full",
##   frequencies = freq[freq > 0 & freq <= pi],
##   ylab = expression(hat(G)[list(n, R)](list(tau[1], tau[2], omega))))


###################################################
### code chunk number 18: quantspec.Rnw:1025-1026 (eval = FALSE)
###################################################
## help("plot-SmoothedPG")


###################################################
### code chunk number 19: quantspec.Rnw:1083-1084
###################################################
CR <- quantilePG(sp500, levels.1 = c(0.05, 0.5, 0.95),
  type = "clipped", type.boot = "mbb", B = 250, l = 32)
freq <- getFrequencies(CR)
plot(CR, levels = c(0.05, 0.5, 0.95),
  frequencies = freq[freq > 0 & freq <= pi],
  ylab = expression({I[list(n, R)]^{list(tau[1], tau[2])}}(omega)))


###################################################
### code chunk number 20: quantspec.Rnw:1093-1094
###################################################
plot(CR, levels = c(0.05, 0.5, 0.95),
  frequencies = freq[freq > 0 & freq <= pi/5],
  ylab = expression({I[list(n, R)]^{list(tau[1], tau[2])}}(omega)))


###################################################
### code chunk number 21: quantspec.Rnw:1107-1108
###################################################
sPG <- smoothedPG(CR, weight = kernelWeight(W = W1, bw = 0.07))
plot(sPG, levels = c(0.05, 0.5, 0.95), type.scaling = "individual",
  frequencies = freq[freq > 0 & freq <= pi], ptw.CIs = 0.1,
  ylab = expression(hat(G)[list(n, R)](list(tau[1], tau[2], omega))))


###################################################
### code chunk number 22: quantspec.Rnw:1118-1119
###################################################
plot(sPG, levels = c(0.05, 0.5, 0.95), type.scaling = "real-imaginary",
  ptw.CIs = 0.1, type.CIs = "boot.full",
  frequencies = freq[freq > 0 & freq <= pi],
  ylab = expression(hat(G)[list(n, R)](list(tau[1], tau[2], omega))))


###################################################
### code chunk number 23: quantspec.Rnw:1158-1159 (eval = FALSE)
###################################################
## help("ts-models")


###################################################
### code chunk number 24: csdQAR (eval = FALSE)
###################################################
## csd <- quantileSD(N = 2^9, seed.init = 2581, type = "copula",
##   ts = ts1, levels.1 = c(0.25, 0.5, 0.75), R = 100, quiet = TRUE)
## plot(csd, ylab = expression(f[list(q[tau[1]], q[tau[2]])](omega)))


###################################################
### code chunk number 25: csdQARhighprec (eval = FALSE)
###################################################
## csd <- quantileSD(N = 2^12, seed.init = 2581, type = "copula",
##   ts = ts1, levels.1 = c(0.25, 0.5, 0.75), R = 50000)
## save(csd, file = "csd-qar1.rdata")


###################################################
### code chunk number 26: quantspec.Rnw:1199-1200 (eval = FALSE)
###################################################
## help("increasePrecision-QuantileSD")


###################################################
### code chunk number 27: csdQARhighprecPlot (eval = FALSE)
###################################################
## load("csd-qar1.rdata")
## plot(csd, frequencies = 2 * pi * (1:2^8) / 2^9,
##   ylab = expression(f[list(q[tau[1]], q[tau[2]])](omega)))


###################################################
### code chunk number 28: quantspec.Rnw:1221-1226 (eval = FALSE)
###################################################
## sCR <- smoothedPG(ts1(512), levels.1 = c(0.25, 0.5, 0.75),
##   weight = kernelWeight(W = W1, bw = 0.1))
## plot(sCR, qsd = csd,
##   ylab = bquote(paste(hat(G)[list(n, R)](list(tau[1], tau[2], omega)),
##   " and ", f[list(q[tau[1]], q[tau[2]])](omega))))


###################################################
### code chunk number 29: quantspec.Rnw:1241-1255 (eval = FALSE)
###################################################
## set.seed(2581)
## ts <- ts1
## N <- 128
## R <- 5000
## 
## freq <- 2  * pi * (1:16) / 32
## levels <- c(0.25, 0.5, 0.75)
## 
## J <- length(freq)
## K <- length(levels)
## 
## sims  <- array(0, dim = c(4, R, J, K, K))
## 
## weight <- kernelWeight(W = W1, bw = 0.3)


###################################################
### code chunk number 30: quantspec.Rnw:1274-1275
###################################################
csd <- quantileSD(N = 2^9, seed.init = 2581, type = "copula",
  ts = ts1, levels.1 = c(0.25, 0.5, 0.75), R = 100, quiet = TRUE)
plot(csd, ylab = expression(f[list(q[tau[1]], q[tau[2]])](omega)))


###################################################
### code chunk number 31: quantspec.Rnw:1300-1306
###################################################
set.seed(020581)
sCR <- smoothedPG(ts1(512), levels.1 = c(0.25, 0.5, 0.75),
  weight = kernelWeight(W = W1, bw = 0.1))
plot(sCR, qsd = csd,
  ylab = bquote(paste(hat(G)[list(n, R)](list(tau[1], tau[2], omega)),
  " and ", f[list(q[tau[1]], q[tau[2]])](omega))))


###################################################
### code chunk number 32: quantspec.Rnw:1319-1330 (eval = FALSE)
###################################################
## for (i in 1:R) {
##   Y <- ts(N)
##   CR  <- quantilePG(Y, levels.1 = levels, type = "clipped")
##   LP  <- quantilePG(Y, levels.1 = levels, type = "qr")
##   sCR <- smoothedPG(CR, weight = weight)
##   sLP <- smoothedPG(LP, weight = weight)
##   sims[1, i, , , ] <- getValues(CR,  frequencies = freq)[, , , 1]
##   sims[2, i, , , ] <- getValues(LP,  frequencies = freq)[, , , 1]
##   sims[3, i, , , ] <- getValues(sCR, frequencies = freq)[, , , 1]
##   sims[4, i, , , ] <- getValues(sLP, frequencies = freq)[, , , 1]
## }


###################################################
### code chunk number 33: quantspec.Rnw:1344-1348 (eval = FALSE)
###################################################
## trueV <- getValues(csd, frequencies = freq)
## SqDev <- array(apply(sims, c(1, 2),
##   function(x) {abs(x - trueV)^2}), dim = c(J, K, K, 4, R))
## rimse <- sqrt(apply(SqDev, c(2, 3, 4), mean))


