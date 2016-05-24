# ex2.02.R
set.seed(123)
theta <- c(0, 5, 3.5)
N <- 100
T <- 1
x <- 10
Z <- rnorm(N)
Dt <- 1/N
A <- theta[3]*sqrt(Dt)*Z
P <- (1-theta[2]*Dt)^(0:(N-1))
X0 <- x
X <- sapply(2:(N+1), function(x) X0*(1-theta[2]*Dt)^(x-1) +
    A[1:(x-1)] %*% P[(x-1):1])
Y <- ts(c(X0,X),start=0, deltat=1/N)
plot(Y)
