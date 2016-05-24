## Active Learning for regression by Integrated Expected
## Conditional Improvement (IECI) with a *un*known constraint
## via Particle Learning on a simple random 1-d function

## load the plgp library
library(plgp)
library(tgp)
library(akima)
library(ellipse)
library(splancs)

## close down old graphics windows and clear session
graphics.off()
rm(list=ls())

## for calculating 2-d data under a unknown (but not hidden) constraint
f2dc <- function(x, y=NULL, epoly)
  {
    ## deal with vector input
    if(is.null(y)) {
      if(!is.matrix(x)) x <- matrix(x, ncol=2)
      y <- x[,2]; x <- x[,1]
    }

    ## first calculate the response
    g <- function(z) {
      return(exp(-(z-1)^2) + exp(-0.8*(z+1)^2) - 0.05*sin(8*(z+0.1)))
    }
    z <- -g(x)*g(y)

    ## then check for being inside the polygon (ellipse)
    cl <- rep(1, length(z))
    if(!is.null(epoly)) {
      xy <- as.points(x,y)
      io <- inout(xy, epoly)
      cl[!io] <- 0
    }

    ## done
    return(data.frame(y=z, c=cl+1))
  }

## describing the constraint region with an ellipse
epoly <- ellipse(rbind(c(1,-0.5),c(-0.5,1)), scale=c(0.75,0.75))
formals(f2dc)$epoly <- epoly

## set bounding rectangle for adaptive sampling
rect <-  rbind(c(-2,2),c(-2,2))
formals(data.ConstGP.improv)$rect <- rect

## set up the generation (Y,C)s for Xs
formals(data.ConstGP.improv)$f <- f2dc

## default ConstGP prior with isotropic covariances
prior <- prior.ConstGP(2, "isotropic")
formals(data.ConstGP.improv)$prior <- prior

## set up the IECI adaptive sampling protocol
formals(data.ConstGP.improv)$cands <- 100

## use akima for interp
library(akima)
formals(data.ConstGP.improv)$interp <- interp

## set up the start and end times
start <- 25
end <- 125

## run PL
out <- PL(dstream=data.ConstGP.improv, ## adaptive design PL via IECI
          start=start, end=end,
          init=draw.ConstGP,  ## init by Metropolis-Hastings
          lpredprob=lpredprob.ConstGP,
          propagate=propagate.ConstGP, prior=prior,
          addpall=addpall.ConstGP, params=params.ConstGP)

## design a grid of predictive locations
XX <- dopt.gp(400, Xcand=lhs(400*10, rect))$XX
XXs <- rectscale(XX, rect)

## sample from the particle posterior predictive distribution 
outp <- papply(XX=XXs, fun=pred.ConstGP, prior=prior)

## extract the mean of the particles predictive
m <- rep(0, nrow(as.matrix(XX)))
for(i in 1:length(outp)) m <- m + outp[[i]]$m
m <- m / length(outp)

## put X back on original scale
X <- rectunscale(PL.env$pall$X, rect)

## set up to make two plots
par(mfrow=c(1,2))

## plot the summary stats of the predictive distribution
## pdf("nanprob_avg.pdf", width=5.5, height=5.5)
image(interp.loess(XX[,1], XX[,2], m),
      main="pred mean", xlab="x1", ylab="x2")
lines(epoly, type="l",lwd=2, lty=2)
points(X)

## mean of the probability of the class 1 label
c1 <- rep(0, nrow(as.matrix(XX)))
for(i in 1:length(outp)) c1 <- c1 + outp[[i]]$class.1
c1 <- c1 / length(outp)

## plot average class.1 probability
## pdf("nanprob_constavg.pdf", width=5.5, height=5.5)
image(interp.loess(XX[,1], XX[,2], c1),
      main="probability of constraint violation",
      xlab="x1", ylab="x2")
lines(epoly, type="l",lwd=2, lty=2)
points(X)

## plot the progress
dev.new()
par(mfrow=c(1,2)) ## two plots
## pdf("nanprob_samps.pdf", width=5.5, height=5.5)
## plot the samples
plot(X[,1], xlab="t", ylim=c(-2,2), main="samples & oracles",
     ylab="x1 & x2")
abline(v=start+1, lwd=2, col=2, lty=2)
points(X[,2], col="blue", xlab="t")
points((start+1):end, PL.env$psave$xstar[,1], col=5, pch=18)
points((start+1):end, PL.env$psave$xstar[,2], col=6, pch=18)
## plot the progres meter
## pdf("nanprob_progress.pdf", width=5.5, height=5.5)
plot(1:end, c(rep(NA, start), PL.env$psave$max.as), type="l", lwd=2,
     xlab="t", ylab="E(reduction in improv)", main="progress meter")
abline(v=start+1, lwd=2, col=2, lty=2)
