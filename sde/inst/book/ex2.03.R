# ex2.03.R
OU <- function(a,b, x, N=1000){
  Y <- numeric(N+1)
  Y[1] <- x
  Z <- rnorm(N)
  Dt <- 1/N
  for(i in 1:N)
    Y[i+1] <- Y[i] - a*Y[i]*Dt + b*sqrt(Dt)*Z[i]
   invisible(Y)
}

OU.vec <- function(a, b, x, N=1000){
  Dt <- 1/N
  Z <- rnorm(N)
  A <- b*sqrt(Dt)*Z
  P <- (1-a*Dt)^(0:(N-1))
   X0 <- x
   X <- c(X0, sapply(2:(N+1), 
      function(x) X0*(1-a*Dt)^(x-1) +
      sum(A[1:(x-1)] * P[(x-1):1])))
 invisible(X)
}

set.seed(123)
system.time(OU(10,5,3.5))
set.seed(123)
system.time(OU.vec(10,5,3.5))
