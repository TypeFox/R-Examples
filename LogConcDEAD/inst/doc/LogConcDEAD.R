### R code from vignette source 'LogConcDEAD.Rnw'

###################################################
### code chunk number 1: LogConcDEAD.Rnw:127-132
###################################################
options(prompt="R> ")
require( "LogConcDEAD" )
require( "logcondens" )
require( "mvtnorm" )
options(width=72)


###################################################
### code chunk number 2: LogConcDEAD.Rnw:417-419 (eval = FALSE)
###################################################
## install.packages("logcondens")
## library("logcondens")


###################################################
### code chunk number 3: setn
###################################################
n <- 100


###################################################
### code chunk number 4: LogConcDEAD.Rnw:430-435
###################################################
set.seed(1)
n <- 100
x <- sort(rgamma(n,shape=2))
out1 <- activeSetLogCon(x) ## logcondens estimate
out2 <- mlelcd(x) ## LogConcDEAD estimate


###################################################
### code chunk number 5: plot:1d (eval = FALSE)
###################################################
## ylim <- c(0,0.4) 
## lgdtxt <- c("LogConcDEAD", "logcondens", "true")
## lgdlty <- c(1,2,3)
## plot(out2, ylim=ylim,lty=1)
## lines(x, exp(out1$phi), lty=2)
## lines(x, x*exp(-x), lty=3)
## legend(x=3, y=0.4, lgdtxt, lty=lgdlty)


###################################################
### code chunk number 6: plot:log1d (eval = FALSE)
###################################################
## ylim <- c(-4,-1)
## lgdtxt <- c("LogConcDEAD", "logcondens", "true")
## lgdlty <- c(1,2,3)
## plot(out2, uselog=TRUE, lty=1)
## lines(x, out1$phi, lty=2)
## lines(x, log(x)-x, lty=3)
## legend(x=3, y=-1, lgdtxt, lty=lgdlty) 


###################################################
### code chunk number 7: fig_1d
###################################################
ylim <- c(0,0.4) 
lgdtxt <- c("LogConcDEAD", "logcondens", "true")
lgdlty <- c(1,2,3)
plot(out2, ylim=ylim,lty=1)
lines(x, exp(out1$phi), lty=2)
lines(x, x*exp(-x), lty=3)
legend(x=3, y=0.4, lgdtxt, lty=lgdlty)


###################################################
### code chunk number 8: fig_log1d
###################################################
ylim <- c(-4,-1)
lgdtxt <- c("LogConcDEAD", "logcondens", "true")
lgdlty <- c(1,2,3)
plot(out2, uselog=TRUE, lty=1)
lines(x, out1$phi, lty=2)
lines(x, log(x)-x, lty=3)
legend(x=3, y=-1, lgdtxt, lty=lgdlty) 


###################################################
### code chunk number 9: setn2
###################################################
n <- 100


###################################################
### code chunk number 10: LogConcDEAD.Rnw:503-507
###################################################
set.seed(22)
d <- 2
n <- 100
x <- matrix(rnorm(n*d),ncol=d)


###################################################
### code chunk number 11: LogConcDEAD.Rnw:519-520
###################################################
out <- mlelcd(x,verbose=50)


###################################################
### code chunk number 12: LogConcDEAD.Rnw:566-567
###################################################
names(out)


###################################################
### code chunk number 13: LogConcDEAD.Rnw:575-576
###################################################
out$logMLE[1:5]


###################################################
### code chunk number 14: LogConcDEAD.Rnw:588-589
###################################################
out$triang[1:5,]


###################################################
### code chunk number 15: LogConcDEAD.Rnw:596-598
###################################################
out$b[1:5,]
out$beta[1:5]


###################################################
### code chunk number 16: LogConcDEAD.Rnw:617-619
###################################################
out$NumberOfEvaluations
out$MinSigma


###################################################
### code chunk number 17: LogConcDEAD.Rnw:623-624
###################################################
out$chull[1:5,]


###################################################
### code chunk number 18: LogConcDEAD.Rnw:630-632
###################################################
out$outnorm[1:5,]
out$outoffset[1:5,]


###################################################
### code chunk number 19: LogConcDEAD.Rnw:662-665
###################################################
g <- interplcd(out, gridlen=200)
g1 <- interpmarglcd(out, marg=1)
g2 <- interpmarglcd(out, marg=2)


###################################################
### code chunk number 20: plot:2d (eval = FALSE)
###################################################
## par(mfrow=c(1,2), pty="s", cex=0.7) #square plots
## plot(out,g=g,addp=FALSE,asp=1)
## plot(out,g=g,uselog=TRUE,addp=FALSE,asp=1)


###################################################
### code chunk number 21: LogConcDEAD.Rnw:679-682
###################################################
png(file="LogConcDEAD-fig_2d.png")
par(mfrow=c(1,2), pty="s", cex=0.7) #square plots
plot(out,g=g,addp=FALSE,asp=1)
plot(out,g=g,uselog=TRUE,addp=FALSE,asp=1)
dev.off()


###################################################
### code chunk number 22: plot:2dmarg (eval = FALSE)
###################################################
## par(mfrow=c(1,2), pty="s", cex=0.7) #normal proportions 
## plot(out,marg=1,g.marg=g1)
## plot(out,marg=2,g.marg=g2)


###################################################
### code chunk number 23: fig_2dmarg
###################################################
par(mfrow=c(1,2), pty="s", cex=0.7) #normal proportions 
plot(out,marg=1,g.marg=g1)
plot(out,marg=2,g.marg=g2)


###################################################
### code chunk number 24: plot:rgl (eval = FALSE)
###################################################
## plot(out,g=g,type="r")


###################################################
### code chunk number 25: plot:rgllog (eval = FALSE)
###################################################
## plot(out,g=g,type="r",uselog=TRUE)


###################################################
### code chunk number 26: LogConcDEAD.Rnw:741-744 (eval = FALSE)
###################################################
## plot(out,g=g,type="r")
## par3d(windowRect = c(55,66,311+256, 322+256))
## rgl.snapshot(file="rglfig.png")


###################################################
### code chunk number 27: LogConcDEAD.Rnw:753-756 (eval = FALSE)
###################################################
## plot(out,g=g,type="r",uselog=TRUE)
## par3d(windowRect = c(55,66,311+256, 322+256))
## rgl.snapshot(file="rgllog.png")


###################################################
### code chunk number 28: LogConcDEAD.Rnw:783-785
###################################################
nsamp <- 1000
mysamp <- rlcd(nsamp,out)


###################################################
### code chunk number 29: LogConcDEAD.Rnw:791-793
###################################################
colMeans(mysamp)
cov(mysamp)


###################################################
### code chunk number 30: LogConcDEAD.Rnw:800-802
###################################################
m <- 10
mypoints <- 1.5*matrix(rnorm(m*d),ncol=d)


###################################################
### code chunk number 31: LogConcDEAD.Rnw:806-807
###################################################
dlcd(mypoints,out)


###################################################
### code chunk number 32: LogConcDEAD.Rnw:822-826
###################################################
myval <- sort(dlcd(mysamp,out))
alpha <- c(.25,.5,.75)
myval[(1-alpha)*nsamp]



###################################################
### code chunk number 33: setn3
###################################################
n <- 100


###################################################
### code chunk number 34: LogConcDEAD.Rnw:841-843 (eval = FALSE)
###################################################
## install.packages("mvtnorm")
## library("mvtnorm")


###################################################
### code chunk number 35: LogConcDEAD.Rnw:846-852
###################################################
set.seed(333)
sigma <- matrix(c(1,0.2,0.2,1),nrow=2)
d <- 2
n <- 100
y <- rmvnorm(n,sigma=0.1*sigma)
xall <- round(y,digits=1)


###################################################
### code chunk number 36: LogConcDEAD.Rnw:865-868
###################################################
tmpw <- getweights(xall)
outw <- mlelcd(tmpw$x,w=tmpw$w)
gw <- interplcd(outw, gridlen=200)


###################################################
### code chunk number 37: plot:bin (eval = FALSE)
###################################################
## par(mfrow=c(1,2), pty="s", cex=0.7) #2 square plots 
## plot(outw,g=gw,asp=1,drawlabels=FALSE)
## plot(outw,g=gw,uselog=TRUE,asp=1,drawlabels=FALSE)


###################################################
### code chunk number 38: LogConcDEAD.Rnw:886-889
###################################################
png(file="LogConcDEAD-fig_bin.png")
par(mfrow=c(1,2), pty="s", cex=0.7) #2 square plots 
plot(outw,g=gw,asp=1,drawlabels=FALSE)
plot(outw,g=gw,uselog=TRUE,asp=1,drawlabels=FALSE)
dev.off()


###################################################
### code chunk number 39: setn4
###################################################
n <- 100


###################################################
### code chunk number 40: LogConcDEAD.Rnw:910-915
###################################################
set.seed(4444)
d <- 3
n <- 100
x <- matrix(rgamma(n*d,shape=2),ncol=d)
out3 <- mlelcd(x)


###################################################
### code chunk number 41: LogConcDEAD.Rnw:922-924
###################################################
mypoints <- c(0,2,4)
dmarglcd(mypoints, out3, marg=1)


###################################################
### code chunk number 42: plot:3deg (eval = FALSE)
###################################################
## par(mfrow=c(2,2),cex=0.8)
## plot(out3, marg=1)
## plot(out3, marg=2)
## plot(out3, marg=3)
## tmp <- seq(min(out3$x), max(out3$x),len=100)
## plot(tmp, dgamma(tmp,shape=2), type="l", 
## xlab="X", ylab="true marginal density")
## title(main="True density")


###################################################
### code chunk number 43: fig_3deg
###################################################
par(mfrow=c(2,2),cex=0.8)
plot(out3, marg=1)
plot(out3, marg=2)
plot(out3, marg=3)
tmp <- seq(min(out3$x), max(out3$x),len=100)
plot(tmp, dgamma(tmp,shape=2), type="l", 
xlab="X", ylab="true marginal density")
title(main="True density")


