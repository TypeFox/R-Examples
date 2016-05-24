data("Somerville", package = "Ecdat")
thetahat <- mean(Somerville$visits)
n <- length(Somerville$visits)
loglnc <- n*thetahat+log(thetahat)*sum(Somerville$visits)-sum(log(factorial(Somerville$visits)))
thetac <- 2.25
loglc <- n*thetac+log(thetac)*sum(Somerville$visits)-sum(log(factorial(Somerville$visits)))
devlnl <- n+sum(Somerville$visits)/thetac

cexdef <- .65
thetahat <- mean(Somerville$visits)
thetac <- 2
loglnc <- -n*thetahat+log(thetahat)*sum(Somerville$visits)-sum(log(factorial(Somerville$visits)))
loglc <- -n*thetac+log(thetac)*sum(Somerville$visits)-sum(log(factorial(Somerville$visits)))
theta <- seq(1.5,3,.05)
logl <- -n*theta+log(theta)*sum(Somerville$visits)-sum(log(factorial(Somerville$visits)))
pdf("graph/threetests.pdf")
plot(logl~theta,type="l",las=1,bty="l",ylim=c(-2950,-2775),xaxs="i",yaxs="i",xlab="",ylab="",cex.axis=cexdef)
lines(c(2.9,2.9),c(-2950,2775))
lines(c(2,2),c(-5989,loglc),lty=3)
lines(c(thetahat,thetahat),c(-5989,loglnc),lty=3)
lines(c(1.5,thetahat),c(loglnc,loglnc),lty=3)
lines(c(1.5,thetac),c(loglc,loglc),lty=3)
devlnl <- -n+sum(Somerville$visits)/thetac
b <- devlnl
a <- loglc-devlnl*2
x0 <- 2
y0 <- a+b*x0
long <- 150
x1 <- x0-sqrt(long/(1+b^2))
x2 <- x0+sqrt(long/(1+b^2))
y1 <- a+b*x1
y2 <- a+b*x2
arrows(x1,y1,x2,y2,code=3,length=.1)
x0 <- 2
y0 <- -2920
x1 <- 1.6
y1 <- -2960
b <- (y1-y0)/(x1-x0)
a <- y0-b*x0
x2 <- 2.6
y2 <- a+b*x2
lines(c(x1,x2),c(y1,y2))
#abline(h=y0,lty=3)
lines(c(2,2.9),rep(a+b*2,2),lty=3)
arrows(thetahat,c(a+b*thetahat),thetahat,-2920,code=3,length=.1)
lines(c(thetahat,2.9),rep(a+b*thetahat,2),lty=3)
text(2.65,-2860,expression(theta-2),cex=cexdef)
text(2.8,-2830,"ln L",cex=cexdef)
text(2.97,a+b*2,0,cex=cexdef)
text(2.97,a+b*thetahat,round(thetahat-2,2),cex=cexdef)
dev.off()
