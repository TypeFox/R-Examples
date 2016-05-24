## Active Learning for regression by Integrated Expected
## Conditional Improvement (IECI) with a *un*known constraint
## via Particle Learning on a simple random 1-d function

## load the plgp library
library(plgp)
library(tgp)

## close down old graphics windows and clear session
graphics.off()
rm(list=ls())

## for calculating 1-d data under a constraint
f1dc <- function(x)
  {
    y <- sin(x) - 2.55*dnorm(x,1.6,0.45)
    y <- y + rnorm(length(x), sd=0.15)
    cl <- as.numeric((x < 2) | (x > 4))+1
    return(data.frame(y=y, c=cl))
  }

## set bounding rectangle for adaptive sampling
rect <- c(0,7)
formals(data.ConstGP.improv)$rect <- rect

## set up the generation Ys for Xs
formals(data.ConstGP.improv)$f <- f1dc

## set the ConstGP prior with tiny nugget for
## small noise data
prior <- prior.ConstGP(2, "isotropic")
prior$grate <- 20
formals(data.ConstGP.improv)$prior <- prior

## set up the IECI adaptive sampling protocol
formals(data.ConstGP.improv)$cands <- 100

## set up the start and end times
start <- 20
end <- 50

## run PL
out <- PL(dstream=data.ConstGP.improv, ## adaptive design PL via IECI
          start=start, end=end,
          init=draw.ConstGP,  ## init with Metropolis-Hastings
          lpredprob=lpredprob.ConstGP, propagate=propagate.ConstGP,
          prior=prior, addpall=addpall.ConstGP,
          params=params.ConstGP)

## design a grid of predictive locations
XX <- seq(rect[1], rect[2], length=101)
XXs <- rectscale(XX, rect)

## sample from the particle posterior predictive distribution 
outp <- papply(XX=XXs, fun=pred.ConstGP, prior=prior)

## put X back on original scale
X <- rectunscale(PL.env$pall$X, rect)

## set up to make four plots
par(mfrow=c(2,2))

## plot the individual lines of the predictive distribution
plot(X, PL.env$pall$Y, main="predictive: each particle", xlab="x", ylab="y")
for(i in 1:length(outp)) {
  lines(XX, outp[[i]]$m)
  lines(XX, outp[[i]]$q1, col=2, lty=2)
  lines(XX, outp[[i]]$q2, col=2, lty=2)
}

## extract the mean of the particles predictive
m <- rep(0, nrow(as.matrix(XX)))
for(i in 1:length(outp)) m <- m + outp[[i]]$m
m <- m / length(outp)

## extract the quantiles of the particles predictive
q2 <- q1 <- rep(0, nrow(as.matrix(XX)))
for(i in 1:length(outp)) {
  q1 <- q1 + outp[[i]]$q1
  q2 <- q2 + outp[[i]]$q2
}
q1 <- q1 / length(outp)
q2 <- q2 / length(outp)

## plot the summary stats of the predictive distribution
## pdf("taddy_avg.pdf", width=5.5, height=5.5)
plot(X, PL.env$pall$Y, main="predictive surface",
     ylim=range(c(q1,q2)), xlab="x", ylab="z")
lines(XX, m, lwd=2)
lines(XX, q1, col=2, lty=2, lwd=2)
lines(XX, q2, col=2, lty=2, lwd=2)
legend("top", c("mean", "90% quantiles"),
       lty=1:2, col=1:2, lwd=2, bty="n")

## plot the individual lines of the posterior classes
plot(XX, outp[[1]]$class.1, main="predictive: Class 1 Prob",
     xlab="x",ylab="prob", xlim=rect, ylim=c(0,1))
for(i in 2:length(outp)) lines(XX, outp[[i]]$class.1)

## mean of the probability of the class 1 label
c1 <- rep(0, nrow(as.matrix(XX)))
for(i in 1:length(outp)) c1 <- c1 + outp[[i]]$class.1
c1 <- c1 / length(outp)

## plot average class.1 probability
## pdf("taddy_constavg.pdf", width=5.5, height=5.5)
plot(XX, c1, type="l", main="probability of constraint violation",
     xlab="x", ylab="prob", lwd=2)

## plot the progress
dev.new()
par(mfrow=c(1,2)) ## two plots
## plot the sampled locations
## pdf("taddy_samps.pdf", width=5.5, height=5.5)
plot(X, xlab="t", ylab="x", main="progress over time")
abline(v=start+1, lwd=2, col=2, lty=2)
abline(h=c(2,4), lwd=2, col=3, lty=3)
points((start+1):end, PL.env$psave$xstar, col=4, pch=18)
## plot the progress meter
## pdf("taddy_progress.pdf", width=5.5, height=5.5)
plot(1:end, c(rep(NA, start), PL.env$psave$max.as), type="l", lwd=2,
     xlab="t", ylab="max log IECI",
     main="progress meter") 
abline(v=start+1, lwd=2, col=2, lty=2)
