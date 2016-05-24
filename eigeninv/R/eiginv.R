eiginv <- function(evals, n, x, y, symmetric=FALSE, stochastic=FALSE) {

soules.1 <- function(x) {
# R code that generates a symmetric Soules matrices
n <- length(x)
tk <- sqrt(cumsum(x * x))

S <- matrix(0, n, n)
S[, 1] <- x / tk[n]
for (j in 2:n) {
	ind <- 1: (n - j + 1)
	S[ind, j] <-  x[ind] * x[n - j + 2] / (tk[n - j + 1] * tk[n - j + 2])
	S[n - j + 2, j] <-  - tk[n - j + 1] / tk[n - j + 2]
}
S
}

soules.2 <- function(x, y=rep(1, length(x))) {
# R code that generates double Soules matrices
# P %*% t(Q) = I
#
n <- length(x)
xy <- sqrt(sum(x * y))
x <- x / xy
y <- y / xy
ss <- 1

P <- Q <- matrix(0, n, n)
P[, 1] <- y 
Q[, 1] <- x 

for (i in 1:(n-1)) {
	cc <- sqrt( 1 - ss * P[n - i + 1, 1] * Q[n - i + 1, 1])
	P[1:(n-i), i + 1] <- ss/cc * Q[n - i + 1, 1]  * P [1:(n-i) ,1] 
	Q[1:(n-i), i + 1] <- ss/cc * P[n - i + 1, 1]  * Q [1:(n-i) ,1] 
	P[n - i + 1, i+1] <- -cc
	Q[n - i + 1, i+1] <- -cc
	ss <- ss / ( 1 - ss * P[n - i + 1, 1] * Q[n - i + 1, 1] )
}
list(P = P, Q = Q)
}


if (missing(evals)) stop("You must provide the eigenvalue spectrum")
if (missing(n)) n <- length(evals)
if (missing(x) & symmetric & stochastic) x <- rep(runif(1), n) 
if (missing(x) & (!stochastic | !symmetric) ) x <- runif(n) 
if (symmetric) {
	ans <- soules.1(x) 
	A <- ans %*% diag(evals) %*% t(ans) 
} else {
	if (missing(y) & stochastic) y <- rep(runif(1), n) 
	if (missing(y) & !stochastic) y <- runif(n) 
	ans <- soules.2(x, y)
	A <- ans$P %*% diag(evals) %*% t(ans$Q)
}
A
}


