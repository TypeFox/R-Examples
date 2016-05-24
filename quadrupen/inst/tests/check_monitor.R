## MONITORING BOUND TO OPTIMALITY
##
## BASICALLY REPRODUCE THE EXAMPLE OF THE PAPER
## GRANDVALET. CHIQUET AND AMBROISE, SPARSITY BY WORST-CASE QUADRATIC PENALTIES
rm(list=ls())
library(quadrupen)

## OTHER TYPE OF DATA
mu <- 3
sigma <- 30
beta <- rep(c(0,-2,2),c(80,10,10))
cor <- 0.75
beta <- rep(c(0,1,0,-1,0), c(25,10,25,10,25))
Soo  <- toeplitz(cor^(0:24)) ## bloc correlation between zero variables
Sww  <- matrix(cor,10,10) ## bloc correlation between active variables
Sigma <- bdiag(Soo,Sww,Soo,Sww,Soo)
p <- length(beta)

diag(Sigma) <- 1
n <- 50
x <- as.matrix(matrix(rnorm(95*n),n,95) %*% chol(Sigma))
y <- 10 + x %*% beta + rnorm(n,0,10)

subset <- sample(1:n, p/2)

lambda2 <- 0.01
fenchel <- elastic.net(x[subset, ], y[subset], lambda2=lambda2, nlambda1=6, control=list(monitor=2))
grandvt <- elastic.net(x[subset, ], y[subset], lambda2=lambda2, nlambda1=6, control=list(monitor=1))

## TRICK FOR FOR LOG-SCALE...
delta_fenchel <- fenchel@monitoring$dist.to.opt
delta_fenchel [delta_fenchel < .Machine$double.eps] <- .Machine$double.eps

delta_grandvt <- grandvt@monitoring$dist.to.opt
delta_grandvt [delta_grandvt < .Machine$double.eps] <- .Machine$double.eps

delta_optimal <- grandvt@monitoring$dist.to.str
delta_optimal[delta_optimal < .Machine$double.eps] <- .Machine$double.eps

plot(log10(delta_optimal), xlab="# iteration", type="l", lty=3, col="black", main="",
     xlim=c(0,length(delta_optimal)), ylim=c(-8, max(log10(c(delta_optimal,delta_fenchel,delta_grandvt)))))
lines(log10(delta_fenchel), xlab="# iteration", lty=2, col="red")
lines(log10(delta_grandvt), xlab="# iteration", lty=1, col="blue")
