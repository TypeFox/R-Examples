## Active Learning for regression by Integrated Expected
## Conditional Improvement (IECI) with a *known* constraint
## via Particle Learning on a simple random 1-d function

## load the plgp library
library(plgp)
library(tgp)

## close down old graphics windows and clear session
graphics.off()
rm(list=ls())

## set up 1-d data
f1d <- function(x, sd=0.15)
  {
    y <- sin(x) - 2.55*dnorm(x,1.6,0.45)
    y <- y + rnorm(length(x), sd=sd)
    return(y)
  }

## set bounding rectangle for adaptive sampling
rect <- c(0,7)
formals(data.GP.improv)$rect <- rect

## set up the generation Ys for Xs
formals(data.GP.improv)$f <- f1d

## reference locations for IECI, thereby encoding
## the *known* constraint
XX <- seq(rect[1], rect[2],length=101)
Xref <- c(XX[XX < 2], XX[XX > 4])
formals(ieci.adapt)$Xref <- Xref

## tiny nugget prior for small noise data

## set the GP prior with tiny nugget for
## small noise data
prior <- prior.GP(1)
prior$grate <- 20
formals(data.GP.improv)$prior <- prior

## set up the IECI adaptive sampling protocol
formals(data.GP.improv)$cands <- 100
formals(data.GP.improv)$adapt <- ieci.adapt

## set up the start and end times
start <- 20
end <- 50

## do the particle learning
out <- PL(dstream=data.GP.improv, ## adaptive design PL via IECI
          start=start, end=end,
          init=draw.GP,  ## init with Metropolis-Hastings
          lpredprob=lpredprob.GP, propagate=propagate.GP,
          addpall=addpall.GP, params=params.GP, prior=prior)

## sample from the particle posterior predictive distribution
XXs <- rectscale(XX, rect)
outp <- papply(XX=XXs, fun=pred.GP, Y=PL.env$pall$Y, quants=TRUE, prior=prior)

## unscale the data locations
X <- rectunscale(PL.env$pall$X, rect)

## set up to make two plots
par(mfrow=c(1,2))

## plot the individual lines of the predictive distribution
plot(X,PL.env$pall$Y, main="predictive: each particle", xlab="x", ylab="y")
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
plot(X, PL.env$pall$Y, main="predictive surface",
     ylim=range(c(q1,q2)), xlab="x", ylab="y")
lines(XX, m, lwd=2)
lines(XX, q1, col=2, lty=2, lwd=2)
lines(XX, q2, col=2, lty=2, lwd=2)
legend("topright", c("mean", "90% quantiles"), lty=1:2, col=1:2, lwd=2)

## plot the progress
dev.new()
par(mfrow=c(1,2)) ## two plots
## plot the sampled locations
plot(X, xlab="t")
abline(v=start+1, lwd=2, col=2, lty=2)
abline(h=c(2,4), lwd=2, col=3, lty=3)
points((start+1):end, PL.env$psave$xstar, col=4, pch=18)
## plot the max log IECI 
plot((start+1):end, PL.env$psave$max.as, type="l", xlab="t",
     ylab="max log IECI", main="progress meter")
