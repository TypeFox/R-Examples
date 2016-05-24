## Active Learning for classification by entropy using
## Particle Learning on a simple 2-d exponential function

## load the plgp library
library(plgp)
library(tgp)
library(akima)

## close down old graphics windows and clear session
graphics.off()
rm(list=ls())

## set bounding rectangle for adaptive sampling
rect <-  rbind(c(-2,2),c(-2,2))
formals(data.CGP.adapt)$rect <- rect

## setup the generation of Cs for Xs
formals(data.CGP.adapt)$f <- exp2d.C

## set candidate locations for the Xs
PL.env$Xcand <- dopt.gp(200, Xcand=lhs(200*10, rect))$XX
formals(data.CGP.adapt)$cands <- NA

## set the CGP prior
prior <- prior.CGP(2, "separable")
formals(data.CGP.adapt)$prior <- prior

## use akima
library(akima)
formals(data.CGP.adapt)$interp <- interp

## start and end
start <- 25
end <- 100

## could be a probem with the first start particles all
## being in the same class

## do the particle learning
out <- PL(dstream=data.CGP.adapt, ## adaptive design PL
          start=start, end=end,
          init=draw.CGP,  ## init with Metropolis-Hastings
          lpredprob.CGP, propagate.CGP, prior=prior,
          addpall.CGP, params.CGP)

## design a grid of predictive locations
XX <- dopt.gp(200, Xcand=lhs(200*10, rect))$XX
XXs <- rectscale(XX, rect)
## truth for comparisons
CC <- exp2d.C(XX)

## sample from the particle posterior predictive distribution
outp <- papply(XX=XXs, fun=pred.CGP, prior=prior)

## extract the mean probability and entropy of the particles predictive
ent <- class <- matrix(NA, nrow=length(outp), ncol=nrow(as.matrix(XX)))
for(i in 1:length(outp)) {
  class[i,] <- apply(outp[[i]], 1, which.max)
  ent[i,] <- apply(outp[[i]], 1, entropy)
}
mclass <- apply(class, 2, mean)
ment <- apply(ent, 2, mean)

## missclassified XXs
CCp <- round(mclass)
miss <- CCp != CC
sum(miss)

## unscale the data locations
X <- rectunscale(PL.env$pall$X, rect)

## plot the summary stats of the predictive distribution
par(mfrow=c(1,2))
## mean
cols <- c(gray(0.85), gray(0.625), gray(0.4))
image(interp(XX[,1], XX[,2], mclass), col=cols,
      xlab="x1", ylab="x2", main="class mean")
points(X); points(XX[miss,], pch=18, col=2)
## entropy
image(interp(XX[,1], XX[,2], ment), 
      xlab="x1", ylab="x2", main="entropy mean")
points(X); points(XX[miss,], pch=18, col=2)
