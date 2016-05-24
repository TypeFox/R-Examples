### R code from vignette source 'dlm.Rnw'

###################################################
### code chunk number 1: dlm.Rnw:26-29
###################################################
options(digits=3)
library(dlm)
set.seed(1963)


###################################################
### code chunk number 2: dlm.Rnw:77-78
###################################################
dlm(FF = 1, V = 0.8, GG = 1, W = 0.1, m0 = 0, C0 = 100)


###################################################
### code chunk number 3: dlm.Rnw:85-86
###################################################
dlmModPoly(order = 1, dV = 0.8, dW = 0.1, C0 = 100)


###################################################
### code chunk number 4: dlm.Rnw:100-101
###################################################
myMod <- dlmModPoly()


###################################################
### code chunk number 5: dlm.Rnw:108-112
###################################################
FF(myMod)
W(myMod)
m0(myMod)
V(myMod) <- 0.8


###################################################
### code chunk number 6: dlm.Rnw:193-194
###################################################
myMod <- dlmModPoly() + dlmModSeas(4)


###################################################
### code chunk number 7: dlm.Rnw:250-253
###################################################
dlmModPoly(dV = 0.2, dW = c(0, 0.5)) %+% 
    (dlmModSeas(4, dV = 0, dW = c(0, 0, 0.35)) + 
     dlmModPoly(1, dV = 0.1, dW = 0.03))


###################################################
### code chunk number 8: dlm.Rnw:286-290
###################################################
u <- rnorm(25)
myMod <- dlmModReg(u, dV = 14.5)
myMod$JFF
head(myMod$X)


###################################################
### code chunk number 9: dlm.Rnw:311-314
###################################################
buildFun <- function(x) {
    dlmModPoly(1, dV = exp(x[1]), dW = exp(x[2]))
}


###################################################
### code chunk number 10: dlm.Rnw:320-325
###################################################
fit <- dlmMLE(Nile, parm = c(0,0), build = buildFun)
fit$conv
dlmNile <- buildFun(fit$par)
V(dlmNile)
W(dlmNile)


###################################################
### code chunk number 11: dlm.Rnw:330-331
###################################################
StructTS(Nile, "level")


###################################################
### code chunk number 12: dlm.Rnw:339-352
###################################################
buildFun <- function(x) {
    m <- dlmModPoly(1, dV = exp(x[1]))
    m$JW <- matrix(1)
    m$X <- matrix(exp(x[2]), nc = 1, nr = length(Nile))
    j <- which(time(Nile) == 1899)
    m$X[j,1] <- m$X[j,1] * (1 + exp(x[3]))
    return(m)
}
fit <- dlmMLE(Nile, parm = c(0,0,0), build = buildFun)
fit$conv
dlmNileJump <- buildFun(fit$par)
V(dlmNileJump)
dlmNileJump$X[c(1, which(time(Nile) == 1899)), 1]


###################################################
### code chunk number 13: figFilter
###################################################
nileJumpFilt <- dlmFilter(Nile, dlmNileJump)
plot(Nile, type = 'o', col = "seagreen")
lines(dropFirst(nileJumpFilt$m), type = 'o', 
      pch = 20, col = "brown")


###################################################
### code chunk number 14: dlm.Rnw:405-413
###################################################
nileJumpFilt <- dlmFilter(Nile, dlmNileJump)
plot(Nile, type = 'o', col = "seagreen")
lines(dropFirst(nileJumpFilt$m), type = 'o', 
      pch = 20, col = "brown")
attach(nileJumpFilt)
v <- unlist(dlmSvd2var(U.C, D.C))
pl <- dropFirst(m) + qnorm(0.05, sd = sqrt(v[-1]))
pu <- dropFirst(m) + qnorm(0.95, sd = sqrt(v[-1]))
detach()
lines(pl, lty = 2, col = "brown") 
lines(pu, lty = 2, col = "brown") 


###################################################
### code chunk number 15: dlm.Rnw:433-440
###################################################
attach(nileJumpFilt)
v <- unlist(dlmSvd2var(U.C, D.C))
pl <- dropFirst(m) + qnorm(0.05, sd = sqrt(v[-1]))
pu <- dropFirst(m) + qnorm(0.95, sd = sqrt(v[-1]))
detach()
lines(pl, lty = 2, col = "brown") 
lines(pu, lty = 2, col = "brown") 


###################################################
### code chunk number 16: figSmooth
###################################################
nileJumpSmooth <- dlmSmooth(nileJumpFilt)
plot(Nile, type = 'o', col = "seagreen")
attach(nileJumpSmooth)
lines(dropFirst(s), type = 'o', pch = 20, col = "brown")
v <- unlist(dlmSvd2var(U.S, D.S))
pl <- dropFirst(s) + qnorm(0.05, sd = sqrt(v[-1]))
pu <- dropFirst(s) + qnorm(0.95, sd = sqrt(v[-1]))
detach()
lines(pl, lty = 2, col = "brown") 
lines(pu, lty = 2, col = "brown") 


###################################################
### code chunk number 17: dlm.Rnw:472-473
###################################################
nileJumpSmooth <- dlmSmooth(nileJumpFilt)
plot(Nile, type = 'o', col = "seagreen")
attach(nileJumpSmooth)
lines(dropFirst(s), type = 'o', pch = 20, col = "brown")
v <- unlist(dlmSvd2var(U.S, D.S))
pl <- dropFirst(s) + qnorm(0.05, sd = sqrt(v[-1]))
pu <- dropFirst(s) + qnorm(0.95, sd = sqrt(v[-1]))
detach()
lines(pl, lty = 2, col = "brown") 
lines(pu, lty = 2, col = "brown") 


###################################################
### code chunk number 18: dlm.Rnw:484-495
###################################################
lGas <- log(UKgas)
dlmGas <- dlmModPoly() + dlmModSeas(4)
buildFun <- function(x) {
    diag(W(dlmGas))[2:3] <- exp(x[1:2])
    V(dlmGas) <- exp(x[3])
    return(dlmGas)
}
(fit <- dlmMLE(lGas, parm = rep(0, 3), build = buildFun))$conv
dlmGas <- buildFun(fit$par)
drop(V(dlmGas))
diag(W(dlmGas))[2:3]


###################################################
### code chunk number 19: dlm.Rnw:502-506
###################################################
gasSmooth <- dlmSmooth(lGas, mod = dlmGas)
x <- cbind(lGas, dropFirst(gasSmooth$s[,c(1,3)]))
colnames(x) <- c("Gas", "Trend", "Seasonal")
plot(x, type = 'o', main = "UK Gas Consumption")


###################################################
### code chunk number 20: dlm.Rnw:511-512
###################################################
plot(x, type = 'o', main = "UK Gas Consumption")


###################################################
### code chunk number 21: dlm.Rnw:524-540
###################################################
gasFilt <- dlmFilter(lGas, mod = dlmGas)
gasFore <- dlmForecast(gasFilt, nAhead = 20)
sqrtR <- sapply(gasFore$R, function(x) sqrt(x[1,1]))
pl <- gasFore$a[,1] + qnorm(0.05, sd = sqrtR)
pu <- gasFore$a[,1] + qnorm(0.95, sd = sqrtR)
x <- ts.union(window(lGas, start = c(1982, 1)), 
              window(gasSmooth$s[,1], start = c(1982, 1)),
              gasFore$a[,1], pl, pu) 
plot(x, plot.type = "single", type = 'o', pch = c(1, 0, 20, 3, 3),
     col = c("darkgrey", "darkgrey", "brown", "yellow", "yellow"),
     ylab = "Log gas consumption")
legend("bottomright", legend = c("Observed", 
                      "Smoothed (deseasonalized)", 
                      "Forecasted level", "90% probability limit"),
       bty = 'n', pch = c(1, 0, 20, 3, 3), lty = 1,
       col = c("darkgrey", "darkgrey", "brown", "yellow", "yellow"))


###################################################
### code chunk number 22: dlm.Rnw:545-553
###################################################
plot(x, plot.type = "single", type = 'o', pch = c(1, 0, 20, 3, 3),
     col = c("darkgrey", "darkgrey", "brown", "yellow", "yellow"),
     ylab = "Log gas consumption")
legend("bottomright", legend = c("Observed", 
                      "Smoothed (deseasonalized)", 
                      "Forecasted level", "90% probability limit"),
       bty = 'n', pch = c(1, 0, 20, 3, 3), lty = 1,
       col = c("darkgrey", "darkgrey", "brown", "yellow", "yellow"))


###################################################
### code chunk number 23: dlm.Rnw:618-622
###################################################
plot(Nile, type = 'o', col = "seagreen")
nileFilt <- dlmFilter(Nile, dlmNile)
for (i in 1:10) # 10 simulated "true" levels 
    lines(dropFirst(dlmBSample(nileFilt)), col = "brown") 


###################################################
### code chunk number 24: dlm.Rnw:626-629
###################################################
plot(Nile, type = 'o', col = "seagreen")
for (i in 1:10) # 10 simulated "true" levels 
    lines(dropFirst(dlmBSample(nileFilt)), col = "brown") 


###################################################
### code chunk number 25: dlm.Rnw:658-662
###################################################
lmixnorm <- function(x, weights, means, sds) {
    log(crossprod(weights, exp(-0.5 * ((x - means) / sds)^2 
                               - log(sds))))
}


###################################################
### code chunk number 26: dlm.Rnw:669-681
###################################################
y <- arms(0, myldens = lmixnorm, 
          indFunc = function(x,...) (x > (-100)) * (x < 100), 
          n = 5000, weights = c(1, 3, 2), 
          means = c(-10, 0, 10), sds = c(7, 5, 2))
summary(y)
library(MASS)
truehist(y, prob = TRUE, ylim = c(0, 0.08), bty = 'o')
curve(colSums(c(1, 3, 2) / 6 *
              dnorm(matrix(x, 3, length(x), TRUE), 
                    mean = c(-10, 0, 10), sd = c(7, 5, 2))), 
      add = TRUE)
legend(-25, 0.07, "True density", lty = 1, bty = 'n')


###################################################
### code chunk number 27: dlm.Rnw:685-690
###################################################
truehist(y, prob = TRUE, ylim = c(0, 0.08), bty = 'o')
curve(colSums(c(1, 3, 2) / 6 *
              dnorm(matrix(x, 3, length(x), TRUE), 
                    mean = c(-10, 0, 10), sd = c(7, 5, 2))), add = TRUE)
legend(-25, 0.07, "True density", lty = 1, bty = 'n')


###################################################
### code chunk number 28: dlm.Rnw:768-773
###################################################
outGibbs <- dlmGibbsDIG(lGas, dlmModPoly(2) + dlmModSeas(4),
                        a.y = 1, b.y = 1000, a.theta = 1, 
                        b.theta = 1000,
                        n.sample = 1100, ind = c(2, 3),
                        save.states = FALSE)


###################################################
### code chunk number 29: dlm.Rnw:778-794
###################################################
burn <- 100
attach(outGibbs)
par(mfrow=c(2,3), mar=c(3.1,2.1,2.1,1.1))
plot(dV[-(1:burn)], type = 'l', xlab="", ylab="", main=expression(sigma^2))
plot(dW[-(1:burn),1], type = 'l', xlab="", ylab="", main=expression(sigma[beta]^2))
plot(dW[-(1:burn),2], type = 'l', xlab="", ylab="", main=expression(sigma[s]^2))
use <- length(dV) - burn
from <- 0.05 * use
at <- pretty(c(0,use),n=3); at <- at[at>=from]
plot(ergMean(dV[-(1:burn)], from), type="l", xaxt="n",xlab="", ylab="")
axis(1, at=at-from, labels=format(at))
plot(ergMean(dW[-(1:burn),1], from), type="l", xaxt="n",xlab="", ylab="")
axis(1, at=at-from, labels=format(at))
plot(ergMean(dW[-(1:burn),2], from), type="l", xaxt="n",xlab="", ylab="")
axis(1, at=at-from, labels=format(at))
detach()


###################################################
### code chunk number 30: dlm.Rnw:805-829
###################################################
burn <- 100
attach(outGibbs)
dV <- dV[-(1:burn)]
dW <- dW[-(1:burn), ]
detach()
par(mfrow=c(2,3), mar=c(3.1,2.1,2.1,1.1))
plot(dV, type = 'l', xlab = "", ylab = "", 
     main = expression(sigma^2))
plot(dW[ , 1], type = 'l', xlab = "", ylab = "", 
     main = expression(sigma[beta]^2))
plot(dW[ , 2], type = 'l', xlab = "", ylab = "", 
     main = expression(sigma[s]^2))
use <- length(dV) - burn
from <- 0.05 * use
at <- pretty(c(0, use), n = 3); at <- at[at >= from]
plot(ergMean(dV, from), type = 'l', xaxt = 'n',
     xlab = "", ylab = "")
axis(1, at = at - from, labels = format(at))
plot(ergMean(dW[ , 1], from), type = 'l', xaxt = 'n',
     xlab = "", ylab = "")
axis(1, at = at - from, labels = format(at))
plot(ergMean(dW[ , 2], from), type = 'l', xaxt = 'n',
     xlab = "", ylab = "")
axis(1, at =  at - from, labels = format(at))


###################################################
### code chunk number 31: dlm.Rnw:836-837
###################################################
mcmcMean(cbind(dV[-(1:burn)], dW[-(1:burn), ]))


###################################################
### code chunk number 32: dlm.Rnw:839-840
###################################################
rm(dV, dW)


