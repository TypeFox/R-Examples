## BOUNDED REGRESSION
##
## CHECK HANDLING OF UNSTABILITY ALONG THE PATH
rm(list=ls())
library(quadrupen)

## Reproduce cases where the quadratic solver fails: high correlation
## for ALL the features (including irrelevant), n = p, too small
## lambda1... (that makes a lot of 'if', but... you never know)
seed <- 8527
cat("\nseed #", seed)
set.seed(seed)
## data generation
mu <- 3
sigma <- 10
beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
Soo  <- matrix(0.75,25,25) ## bloc correlation between zero variables
Sww  <- matrix(0.75,10,10) ## bloc correlation between active variables
Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo) + 0.2
diag(Sigma) <- 1
n <- 100
x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
y <- mu + x %*% beta + rnorm(n,0,10)

## fails of the quadratic solver, correctly handled by 'bullet proof' mode...
fit01 <- bounded.reg(x,y, lambda2=0,min.ratio=1e-4, control=list(bulletproof=FALSE)) ## Early fail in Cholesky decompostion ...
fit11 <- bounded.reg(x,y, lambda2=0,min.ratio=1e-4)
fit02 <- bounded.reg(x,y, lambda2=0.1,min.ratio=1e-4, control=list(bulletproof=FALSE)) ## out of convergence
fit12 <- bounded.reg(x,y, lambda2=0.1,min.ratio=1e-4)
