## Particle Learning for classification Gaussian Processes
## on the Cushings data from Ripley

## load the necessary libraries
library(tgp)
library(akima)
library(MASS)
library(plgp)

## close down old graphics windows and clear session
graphics.off()
rm(list=ls())

## get the data in a bounding rectangle
data(Cushings)
rect <- rbind(c(0.5, 4), c(-3.5, 3))
X <- log(Cushings[1:21,1:2])
C <- as.numeric(Cushings[1:21,3])

## coerse the input matrix and adjust by the rectangle
## and put data into data-generating function
Xs <- rectscale(X, rect)
formals(data.CGP)$X <- Xs
formals(data.CGP)$C <- C

## set up starting and ending times:
## by sepcifying start = end we can do Metropolis-Hastings
## for comparisons
start <- ncol(Xs) + 5*length(unique(C))
end <- nrow(Xs)

## default prior
prior <- prior.CGP(2)

## could be a probem with the first start particles all
## being in the same class

## do the particle learning
out <- PL(dstream=data.CGP, ## static PL
          start=start, end=end,
          init=draw.CGP,  ## initializing with Metropolis-Hastings
          lpredprob.CGP, propagate.CGP, prior=prior,
          addpall.CGP, params.CGP)

## design a grid of predictive locations
ng <- 50
Tetra <- seq(0,4.5,length=ng)
Preg <- seq(-4,3,length=ng)
XX <- expand.grid(Tetra, Preg )
XXs <- rectscale(XX, rect)

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

## unscale the data locations
X <- rectunscale(PL.env$pall$X, rect)
  
## plot the summary stats of the predictive distribution
par(mfrow=c(1,2))
## mean
cols <- c(gray(0.85), gray(0.625), gray(0.4))
image(interp(XX[,1], XX[,2], mclass), col=cols,
      xlab="x1", ylab="x2", main="class mean")
text(X[,1], X[,2], labels=C)
## entropy
image(interp(XX[,1], XX[,2], ment), 
      xlab="x1", ylab="x2", main="entropy mean")
text(X[,1], X[,2], labels=C) 
