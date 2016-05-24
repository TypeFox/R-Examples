### R code from vignette source 'tgp.Rnw'

###################################################
### code chunk number 1: tgp.Rnw:26-28
###################################################
library(tgp)
options(width=65)


###################################################
### code chunk number 2: tgp.Rnw:141-142
###################################################
bgp


###################################################
### code chunk number 3: gpllm
###################################################
hist(c(rgamma(100000,1,20), rgamma(100000,10,10)), 
     breaks=50, xlim=c(0,2), freq=FALSE, ylim=c(0,3),
     main = "p(d) = G(1,20) + G(10,10)", xlab="d")
d <- seq(0,2,length=1000)
lines(d,0.2+0.7/(1+exp(-10*(d-0.5))))
abline(h=1, lty=2)
legend(x=1.25, y=2.5, c("p(b) = 1", "p(b|d)"), lty=c(1,2))


###################################################
### code chunk number 4: tgp.Rnw:624-625
###################################################
graphics.off()


###################################################
### code chunk number 5: linear.iRnw:1-4
###################################################
library(tgp)
##options(width=65)
seed <- 0; set.seed(seed)


###################################################
### code chunk number 6: linear.iRnw:20-24
###################################################
# 1-d linear data input and predictive data
X <- seq(0,1,length=50)  # inputs
XX <- seq(0,1,length=99) # predictive locations
Z <- 1 + 2*X + rnorm(length(X),sd=0.25) # responses


###################################################
### code chunk number 7: linear.iRnw:29-30
###################################################
lin.blm <- blm(X=X, XX=XX, Z=Z)


###################################################
### code chunk number 8: linear-blm
###################################################
plot(lin.blm, main='Linear Model,', layout='surf')
abline(1,2,lty=3,col='blue')


###################################################
### code chunk number 9: linear.iRnw:38-39
###################################################
graphics.off()


###################################################
### code chunk number 10: linear.iRnw:72-73
###################################################
lin.gpllm <- bgpllm(X=X, XX=XX, Z=Z)


###################################################
### code chunk number 11: linear-gplm
###################################################
plot(lin.gpllm, main='GP LLM,', layout='surf')
abline(1,2,lty=4,col='blue')


###################################################
### code chunk number 12: linear.iRnw:81-82
###################################################
graphics.off()


###################################################
### code chunk number 13: linear.iRnw:102-106
###################################################
lin.gpllm.tr <- bgpllm(X=X, XX=0.5, Z=Z, pred.n=FALSE, trace=TRUE,
                       verb=0)
mla <- mean(lin.gpllm.tr$trace$linarea$la)
mla


###################################################
### code chunk number 14: linear.iRnw:112-113
###################################################
1-mean(lin.gpllm.tr$trace$XX[[1]]$b1)


###################################################
### code chunk number 15: sin.iRnw:4-7
###################################################
library(tgp)
##options(width=65)
seed <- 0; set.seed(seed)


###################################################
### code chunk number 16: sin.iRnw:25-31
###################################################
X <- seq(0,20,length=100)
XX <- seq(0,20,length=99)
Ztrue <- (sin(pi*X/5) + 0.2*cos(4*pi*X/5)) * (X <= 9.6)
lin <- X>9.6; 
Ztrue[lin] <- -1 + X[lin]/10
Z <- Ztrue + rnorm(length(Ztrue), sd=0.1)


###################################################
### code chunk number 17: sin.iRnw:36-37
###################################################
sin.bgp <- bgp(X=X, Z=Z, XX=XX, verb=0)


###################################################
### code chunk number 18: sin-bgp
###################################################
plot(sin.bgp, main='GP,', layout='surf')
lines(X, Ztrue, col=4, lty=2, lwd=2)


###################################################
### code chunk number 19: sin.iRnw:45-46
###################################################
graphics.off()


###################################################
### code chunk number 20: sin.iRnw:62-63
###################################################
sin.btlm <- btlm(X=X, Z=Z, XX=XX)


###################################################
### code chunk number 21: sin-btlm
###################################################
plot(sin.btlm, main='treed LM,', layout='surf')
lines(X, Ztrue, col=4, lty=2, lwd=2)


###################################################
### code chunk number 22: sin.iRnw:76-77
###################################################
graphics.off()


###################################################
### code chunk number 23: sin-btlmtrees
###################################################
tgp.trees(sin.btlm)


###################################################
### code chunk number 24: sin.iRnw:84-85
###################################################
graphics.off()


###################################################
### code chunk number 25: sin.iRnw:106-107
###################################################
sin.btgp <- btgp(X=X, Z=Z, XX=XX, verb=0)


###################################################
### code chunk number 26: sin-btgp
###################################################
plot(sin.btgp, main='treed GP,', layout='surf')
lines(X, Ztrue, col=4, lty=2, lwd=2)


###################################################
### code chunk number 27: sin.iRnw:115-116
###################################################
graphics.off()


###################################################
### code chunk number 28: exp.iRnw:4-8
###################################################
library(tgp)
library(maptree)
##options(width=65)
seed <- 0; set.seed(seed)


###################################################
### code chunk number 29: exp.iRnw:26-29
###################################################
exp2d.data <- exp2d.rand()
X <- exp2d.data$X; Z <- exp2d.data$Z
XX <- exp2d.data$XX


###################################################
### code chunk number 30: exp.iRnw:39-40
###################################################
exp.bgp <- bgp(X=X, Z=Z, XX=XX, corr="exp", verb=0) 	


###################################################
### code chunk number 31: exp-bgp
###################################################
plot(exp.bgp, main='GP,')


###################################################
### code chunk number 32: exp.iRnw:47-48
###################################################
graphics.off()


###################################################
### code chunk number 33: exp.iRnw:71-72
###################################################
exp.btgp <- btgp(X=X, Z=Z, XX=XX, corr="exp", verb=0)


###################################################
### code chunk number 34: exp-btgp
###################################################
plot(exp.btgp, main='treed GP,')


###################################################
### code chunk number 35: exp.iRnw:79-80
###################################################
graphics.off()


###################################################
### code chunk number 36: exp-btgptrees
###################################################
tgp.trees(exp.btgp)


###################################################
### code chunk number 37: exp.iRnw:87-88
###################################################
graphics.off()


###################################################
### code chunk number 38: exp.iRnw:112-113
###################################################
exp.btgpllm <- btgpllm(X=X, Z=Z, XX=XX, corr="exp", R=2) 	


###################################################
### code chunk number 39: exp-btgpllm
###################################################
plot(exp.btgpllm, main='treed GP LLM,')


###################################################
### code chunk number 40: exp.iRnw:120-121
###################################################
graphics.off()


###################################################
### code chunk number 41: exp-1dbtgpllm1
###################################################
plot(exp.btgpllm, main='treed GP LLM,', proj=c(1))


###################################################
### code chunk number 42: exp.iRnw:143-144
###################################################
graphics.off()


###################################################
### code chunk number 43: exp-1dbtgpllm2
###################################################
plot(exp.btgpllm, main='treed GP LLM,', proj=c(2))


###################################################
### code chunk number 44: exp.iRnw:150-151
###################################################
graphics.off()


###################################################
### code chunk number 45: moto.iRnw:4-7
###################################################
library(tgp)
##options(width=65)
seed <- 0; set.seed(seed)


###################################################
### code chunk number 46: moto.iRnw:30-33
###################################################
library(MASS)
X <- data.frame(times=mcycle[,1])
Z <- data.frame(accel=mcycle[,2])


###################################################
### code chunk number 47: moto.iRnw:39-40
###################################################
moto.bgp <- bgp(X=X, Z=Z, verb=0)


###################################################
### code chunk number 48: moto-bgp
###################################################
plot(moto.bgp, main='GP,', layout='surf')


###################################################
### code chunk number 49: moto.iRnw:48-49
###################################################
graphics.off()


###################################################
### code chunk number 50: moto.iRnw:62-63
###################################################
moto.btlm <- btlm(X=X, Z=Z, verb=0)


###################################################
### code chunk number 51: moto-btlm
###################################################
plot(moto.btlm, main='Bayesian CART,', layout='surf')


###################################################
### code chunk number 52: moto.iRnw:72-73
###################################################
graphics.off()


###################################################
### code chunk number 53: moto.iRnw:90-92
###################################################
moto.btgpllm <- btgpllm(X=X, Z=Z, bprior="b0", verb=0)
moto.btgpllm.p <- predict(moto.btgpllm) ## using MAP


###################################################
### code chunk number 54: moto-btgp
###################################################
par(mfrow=c(1,2))
plot(moto.btgpllm, main='treed GP LLM,', layout='surf')
plot(moto.btgpllm.p, center='km', layout='surf')


###################################################
### code chunk number 55: moto.iRnw:103-104
###################################################
graphics.off()


###################################################
### code chunk number 56: moto-btgpq
###################################################
par(mfrow=c(1,2))
plot(moto.btgpllm, main='treed GP LLM,', layout='as')
plot(moto.btgpllm.p, as='ks2', layout='as')


###################################################
### code chunk number 57: moto.iRnw:114-115
###################################################
graphics.off()


###################################################
### code chunk number 58: fried.iRnw:4-7
###################################################
library(tgp)
##options(width=65)
seed <- 0; set.seed(seed)


###################################################
### code chunk number 59: fried.iRnw:38-42
###################################################
f <- friedman.1.data(200)
ff <- friedman.1.data(1000)
X <- f[,1:10]; Z <- f$Y
XX <- ff[,1:10]


###################################################
### code chunk number 60: fried.iRnw:49-52
###################################################
fr.btlm <- btlm(X=X, Z=Z, XX=XX, tree=c(0.95,2), pred.n=FALSE, verb=0)
fr.btlm.mse <- sqrt(mean((fr.btlm$ZZ.mean - ff$Ytrue)^2))
fr.btlm.mse


###################################################
### code chunk number 61: fried.iRnw:55-58
###################################################
fr.bgpllm <- bgpllm(X=X, Z=Z, XX=XX, pred.n=FALSE, verb=0)
fr.bgpllm.mse <- sqrt(mean((fr.bgpllm$ZZ.mean - ff$Ytrue)^2))
fr.bgpllm.mse


###################################################
### code chunk number 62: fried.iRnw:67-70
###################################################
XX1 <- matrix(rep(0,10), nrow=1)
fr.bgpllm.tr <- bgpllm(X=X, Z=Z, XX=XX1, pred.n=FALSE, trace=TRUE, 
                       m0r1=FALSE, verb=0)


###################################################
### code chunk number 63: fried.iRnw:80-82
###################################################
trace <- fr.bgpllm.tr$trace$XX[[1]]
apply(trace[,27:36], 2, mean)


###################################################
### code chunk number 64: fried.iRnw:88-89
###################################################
mean(fr.bgpllm.tr$trace$linarea$ba)


###################################################
### code chunk number 65: fried.iRnw:95-96
###################################################
summary(trace[,9:10])


###################################################
### code chunk number 66: fried.iRnw:99-100
###################################################
apply(trace[,11:15], 2, mean)


###################################################
### code chunk number 67: as.iRnw:4-8
###################################################
library(tgp)
library(maptree)
#options(width=65)
seed <- 0; set.seed(seed)


###################################################
### code chunk number 68: as.iRnw:14-18
###################################################
exp2d.data <- exp2d.rand(lh=0, dopt=10)
X <- exp2d.data$X
Z <- exp2d.data$Z
Xcand <- lhs(1000, rbind(c(-2,6),c(-2,6)))


###################################################
### code chunk number 69: as.iRnw:31-32
###################################################
exp1 <- btgpllm(X=X, Z=Z, pred.n=FALSE, corr="exp", R=5, verb=0)


###################################################
### code chunk number 70: as-mapt
###################################################
tgp.trees(exp1)


###################################################
### code chunk number 71: as.iRnw:39-40
###################################################
graphics.off()


###################################################
### code chunk number 72: as.iRnw:55-57
###################################################
XX <- tgp.design(200, Xcand, exp1)
XX <- rbind(XX, c(-sqrt(1/2),0))


###################################################
### code chunk number 73: as-cands
###################################################
plot(exp1$X, pch=19, cex=0.5)
points(XX)
mapT(exp1, add=TRUE)


###################################################
### code chunk number 74: as.iRnw:74-75
###################################################
graphics.off()


###################################################
### code chunk number 75: as.iRnw:89-91
###################################################
exp.as <- btgpllm(X=X, Z=Z, XX=XX, corr="exp", improv=TRUE, 
                        Ds2x=TRUE, R=5, verb=0)


###################################################
### code chunk number 76: as-expas
###################################################
par(mfrow=c(1,3), bty="n")
plot(exp.as, main="tgpllm,", layout="as", as="alm")
plot(exp.as, main="tgpllm,", layout='as', as='alc')
plot(exp.as, main="tgpllm,", layout='as', as='improv')


###################################################
### code chunk number 77: as.iRnw:109-110
###################################################
graphics.off()


###################################################
### code chunk number 78: traces.iRnw:4-7
###################################################
library(tgp)
##options(width=65)
seed <- 0; set.seed(seed)


###################################################
### code chunk number 79: traces.iRnw:43-47
###################################################
exp2d.data <- exp2d.rand(n2=150, lh=0, dopt=10)
X <- exp2d.data$X
Z <- exp2d.data$Z
XX <- rbind(c(0,0),c(2,2),c(4,4))


###################################################
### code chunk number 80: traces.iRnw:52-55
###################################################
out <- btgpllm(X=X, Z=Z, XX=XX, corr="exp", bprior="b0", 
               pred.n=FALSE, Ds2x=TRUE, R=10, 
               trace=TRUE, verb=0)


###################################################
### code chunk number 81: traces.iRnw:59-60
###################################################
out$trace


###################################################
### code chunk number 82: traces-XXd
###################################################
trXX <- out$trace$XX; ltrXX <- length(trXX)
y <- trXX[[1]]$d
for(i in 2:ltrXX) y <- c(y, trXX[[i]]$d)
plot(log(trXX[[1]]$d), type="l", ylim=range(log(y)), ylab="log(d)",
     main="range (d) parameter traces")
names <- "XX[1,]"
for(i in 2:ltrXX) {
  lines(log(trXX[[i]]$d), col=i, lty=i)
  names <- c(names, paste("XX[", i, ",]", sep=""))
}
legend("bottomleft", names, col=1:ltrXX, lty=1:ltrXX)


###################################################
### code chunk number 83: traces.iRnw:92-93
###################################################
graphics.off()


###################################################
### code chunk number 84: traces.iRnw:110-112
###################################################
linarea <- mean(out$trace$linarea$la)
linarea


###################################################
### code chunk number 85: traces-la
###################################################
hist(out$trace$linarea$la)


###################################################
### code chunk number 86: traces.iRnw:119-120
###################################################
graphics.off()


###################################################
### code chunk number 87: traces.iRnw:135-141
###################################################
m <- matrix(0, nrow=length(trXX), ncol=3)#ncol=5)
for(i in 1:length(trXX))
  m[i,] <- as.double(c(out$XX[i,], mean(trXX[[i]]$b)))
m <- data.frame(cbind(m, 1-m[,3]))
names(m)=c("XX1","XX2","b","pllm")
m


###################################################
### code chunk number 88: traces-alc
###################################################
trALC <- out$trace$preds$Ds2x
y <- trALC[,1]
for(i in 2:ncol(trALC)) y <- c(y, trALC[,i])
plot(log(trALC[,1]), type="l", ylim=range(log(y)), ylab="Ds2x",
     main="ALC: samples from Ds2x")
names <- "XX[1,]"
for(i in 2:ncol(trALC)) {
  lines(log(trALC[,i]), col=i, lty=i)
  names <- c(names, paste("XX[", i, ",]", sep=""))
}
legend("bottomright", names, col=1:ltrXX, lty=1:ltrXX)


###################################################
### code chunk number 89: traces.iRnw:162-163
###################################################
graphics.off()


###################################################
### code chunk number 90: pred.iRnw:4-8
###################################################
library(tgp)
library(maptree)
##options(width=65)
seed <- 0; set.seed(seed)


###################################################
### code chunk number 91: pred.iRnw:18-23
###################################################
library(MASS)
out <- btgpllm(X=mcycle[,1], Z=mcycle[,2], bprior="b0", 
	       pred.n=FALSE, verb=0) 
save(out, file="out.Rsave")
out <- NULL


###################################################
### code chunk number 92: pred.iRnw:32-35
###################################################
load("out.Rsave")
XX <- seq(2.4, 56.7, length=200)
out.kp <- predict(out, XX=XX, pred.n=FALSE)


###################################################
### code chunk number 93: pred.iRnw:40-41
###################################################
out.p <- predict(out, XX=XX, pred.n=FALSE, BTE=c(0,1000,1))


###################################################
### code chunk number 94: pred.iRnw:50-51
###################################################
out2 <- predict(out, XX, pred.n=FALSE, BTE=c(0,2000,2), MAP=FALSE)


###################################################
### code chunk number 95: pred-kp
###################################################
plot(out.kp, center="km", as="ks2")


###################################################
### code chunk number 96: pred.iRnw:62-63
###################################################
graphics.off()


###################################################
### code chunk number 97: pred-p
###################################################
plot(out.p)


###################################################
### code chunk number 98: pred.iRnw:70-71
###################################################
graphics.off()


###################################################
### code chunk number 99: pred-2
###################################################
plot(out2)


###################################################
### code chunk number 100: pred.iRnw:78-79
###################################################
graphics.off()


###################################################
### code chunk number 101: pred.iRnw:102-103
###################################################
unlink("out.Rsave")


