###############################################################################
## Optimally Robust Estimation
###############################################################################

###############################################################################
## Example 1
###############################################################################
require(MASS)
require(robustbase)
library(ISwR)
data(thuesen)
attach(thuesen)

## LS-estimator
fit.LS <- lm(short.velocity ~ blood.glucose)
## M-estimator
fit.M <- rlm(short.velocity ~ blood.glucose)
## MM-estimator
fit.MM <- lmrob(short.velocity ~ blood.glucose)
## LTS-estimator
fit.LTS <- ltsReg(short.velocity ~ blood.glucose)

## AL-estimators
## regressor distribution: design measure
require(RobRex)
K <- DiscreteMVDistribution(cbind(1,blood.glucose))
## ALc-estimator: conditional neighborhood; i.e., only y is contaminated
## about 10 sec. on my system
system.time(IC1 <- rgsOptIC.ALc(r = 0.5, K = K, theta = fit.M$coeff, scale = fit.M$s))
ALc1 <- oneStepEstimator(cbind(1,as.matrix(thuesen)), IC1, c(fit.M$coeff, fit.M$s))
## AL-estimator: unconditional neighborhood; i.e., x and y are contaminated
## about 85 sec. on my system
system.time(IC2 <- rgsOptIC.AL(r = 0.5, K = K, theta = fit.MM$coeff, scale = fit.MM$s))
AL1 <- oneStepEstimator(cbind(1,as.matrix(thuesen)), IC2, c(fit.MM$coeff, fit.MM$s))

## Plot
require(RColorBrewer)
myCol <- brewer.pal(6, "Set1")
plot(short.velocity ~ blood.glucose, ylab = "fasting blood glucose [mmol/l]", 
     xlab = "mean circumferential shortening velocity [%/s]", 
     main = "Ventricular shortening velocity", pch = 20)
abline(fit.LS, lwd = 2, col = myCol[1])
abline(fit.M, lwd = 2, col = myCol[2])
abline(fit.MM, lwd = 2, col = myCol[3])
abline(fit.LTS, lwd = 2, col = myCol[4])
lines(c(1, c(blood.glucose,25)), ALc1[1] + ALc1[2]*c(1,c(blood.glucose,25)), 
      col = myCol[5], lwd = 2)
lines(c(1, c(blood.glucose,25)), AL1[1] + AL1[2]*c(1,c(blood.glucose,25)), 
      col = myCol[6], lwd = 1)
legend("topleft", legend = c("LS", "M", "MM", "LTS", "ALc", "AL"), 
       fill = myCol, ncol = 2)
detach(thuesen)

###############################################################################
## Example 2
###############################################################################
data(phones)
attach(phones)

## LS estimator
fit2.LS <- lm(calls ~ year)
## M estimator
fit2.M <- rlm(calls ~ year, maxit = 50)
## MM estimator
fit2.MM <- lmrob(calls ~ year)
## LTS estimator
fit2.LTS <- ltsReg(calls ~ year)

## AL estimators
## regressor distribution: design measure
K <- DiscreteMVDistribution(cbind(1,year))
## ALc estimator: only y contaminated
system.time(IC1 <- rgsOptIC.ALc(r = 0.5, K = K, theta = fit2.M$coeff, scale = fit2.M$s))
## takes about 9 sec. on my system
ALc2 <- oneStepEstimator(cbind(1,year,calls), IC1, c(fit2.M$coeff, fit2.M$s))
## AL estimator: x and y contaminated
system.time(IC2 <- rgsOptIC.AL(r = 0.5, K = K, theta = fit2.MM$coeff, scale = fit2.MM$s))
## takes about 80 sec. on my system
AL2 <- oneStepEstimator(cbind(1,year,calls), IC2, c(fit2.MM$coeff, fit2.MM$s))

## Plot
plot(calls ~ year, ylab = "phone calls [Mio.]", xlab = "year", 
     main = "Belgium Phone Calls 1950-1973", pch = 20)
abline(fit2.LS, lwd = 2, col = "black", col = myCol[1])
abline(fit2.M, lwd = 2, col = myCol[2])
abline(fit2.MM, lwd = 2, col = myCol[3])
abline(fit2.LTS, lwd = 2, col = myCol[4], lty = 2)
lines(c(1, c(year,75)), ALc2[1] + ALc2[2]*c(1,c(year,75)), lwd = 2, col = myCol[5])
lines(c(1, c(year,75)), AL2[1] + AL2[2]*c(1,c(year,75)), lwd = 2, col = myCol[6])
legend("topleft", legend = c("LS", "M", "MM", "LTS", "ALc", "AL"), 
       fill = myCol, ncol = 2)
