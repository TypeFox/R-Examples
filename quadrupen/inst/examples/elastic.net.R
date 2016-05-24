rm(list=ls())
library(quadrupen)
## Simulating multivariate Gaussian with blockwise correlation
## and piecewise constant vector of parameters
beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
Soo  <- matrix(0.75,25,25) ## bloc correlation between zero variables
Sww  <- matrix(0.75,10,10) ## bloc correlation between active variables
Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo) + 0.2
diag(Sigma) <- 1
n <- 100
x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
y <- 10 + x %*% beta + rnorm(n,0,10)

## Structured Elastic.net
p <- ncol(x)
C <- bandSparse(p,k=0:1,diagonals=list(rep(1,p),rep(-1,p-1)))

## This gives a great advantage to the elastic-net
## for support recovery
beta.lasso <- slot(crossval(x,y,lambda2=0) , "beta.min")
beta.enet  <- slot(crossval(x,y,lambda2=10), "beta.min")
beta.struc <- slot(crossval(x,y,lambda2=10,struct=t(C) %*% C), "beta.min")

cat("\nFalse positives for the Lasso:", sum(sign(beta) != sign(beta.lasso)))
cat("\nFalse positives for the Elastic-net:", sum(sign(beta) != sign(beta.enet)))
cat("\nFalse positives for Structured Elastic-net:", sum(sign(beta) != sign(beta.struc)))
cat("\nDONE.\n")

## Comparing the solution path pf the LASSO, the Elastic-net and the
## Structured Elastic-net
labels <- rep("irrelevant", length(beta))
labels[beta != 0] <- "relevant"
plot(elastic.net(x,y,lambda2=0), label=labels) ## a mess
plot(elastic.net(x,y,lambda2=10), label=labels) ## a lot better
plot(elastic.net(x,y,lambda2=10,struct=solve(Sigma)), label=labels) ## even better
