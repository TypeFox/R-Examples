## CROSS-VALIDATION OF BOUNDED REGRESSION
##
## CHECK CORRECT BEHAVIOR WHEN THE UNDERLYING BOUNDED.REG FUNCTION
## EXPERIENCES UNSTABILITY

library(quadrupen)
set.seed(111)
## mutivariate Gaussian data
mu <- 3
sigma <- 10
n <- 100
beta <- rep(c(0,-2,2),c(80,10,10))
cor <- 0.8
eps <- 0.25 # correlation between relevant and irrelevant variable
S11 <- toeplitz(cor^(0:(80-1))) ## Toeplitz correlation for irrelevant variables
S22 <- matrix(cor,10,10) ## bloc correlation bewteen relevant variables
diag(S22) <- 1
Sigma <- bdiag(S11,S22,S22) + eps

x <- as.matrix(matrix(rnorm(100*n),n,100) %*% chol(Sigma))
y <- mu + x %*% beta + rnorm(n, 0, sigma)

## THESE SETTINGS SHOULD INDUCE EARLY STOPS: CHECK RELEVANCE OF THE
## CROSSVAL FUNCTION

## FIRST, with bulletproof mode: should go at the end of the path, but slow (proximal)
cv.simple.bp  <- crossval(x, y, penalty="bounded.reg", lambda2=0, min.ratio=1e-5)
plot(cv.simple.bp)

## SECOND, without bulletproof: early stop, fast, but possible larger
## standard errors because of less points to average for small lambdas (irrelevant anyway in a sparse setting!)
cv.simple.nbp <- crossval(x, y, penalty="bounded.reg", lambda2=0, min.ratio=1e-5, control=list(bulletproof=FALSE))
plot(cv.simple.nbp)

