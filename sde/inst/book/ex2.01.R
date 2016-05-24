# ex2.01.R
set.seed(123)
N <- 100
T <- 1
x <- 10
theta <- c(0, 5, 3.5)
Dt <- 1/N
Y <- numeric(N+1)
Y[1] <- x
Z <- rnorm(N)
for(i in 1:N)
  Y[i+1] <- Y[i] + (theta[1] - theta[2]*Y[i])*Dt + theta[3]*sqrt(Dt)*Z[i]
Y <- ts(Y,start=0, deltat=1/N)
plot(Y)
