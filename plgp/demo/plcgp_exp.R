## Classification Gaussian Processes
## on a simple 2-d exponential function

## load the necessary libraries
library(plgp)
library(tgp)
library(akima)

## close down old graphics windows and clear session
graphics.off()
rm(list=ls())

## create the design and data in a bounding rectangle
rect <-  rbind(c(-2,2),c(-2,2))
X <- dopt.gp(125, Xcand=lhs(10*125, rect))$XX
C <- exp2d.C(X)

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
## being in the same class -- causes an error

## Particle Learning Inference!
out <- PL(dstream=data.CGP, ## static PL
          start=start, end=end,
          init=draw.CGP,  ## initializing with Metropolis-Hastings
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

## look at a historgram of the parameters
params <- params.CGP()
dev.new()
par(mfrow=c(3,2))
hist(params$d.2); hist(params$d.3)
hist(params$g.2); hist(params$g.3)
hist(params$lpost.2); hist(params$lpost.3)
