Q <- matrix(c(-24, 2, 0, 24, -14, 4, 0, 12, -4), 3, 3)
library(Matrix)
(Phalf <- expm(Q/2))
(Pone <- expm(Q))
(Ptwo <- expm(2*Q))

Phalf %*% Phalf
Pone %*% Pone

pi0 <- c(0,0,1)
(pi0 %*% Phalf)
(pi0 %*% Pone)
(pi0 %*% Ptwo)

A <- t(Q)
A[1,] <- c(1,1,1)
solve(A, c(1,0,0))

# exp time to hit 1
TQ <- Q[-1,-1]
solve(TQ, rep(-1,2))

# via simulation
(rates <- -diag(Q))
(P <- diag(1/rates) %*% Q + diag(1,3))

hit <- function(rates, P) {
  x <- 3
  t <- 0
  while (x > 1) {
    t <- t + rexp(1, rates[x])
    x <- sample(1:3, 1, prob = P[x,])
    #cat(x, t, "\n")
  }
  return(t)
}

N <- 1000000
hits <- rep(0, N)
for (i in 1:N) hits[i] <- hit(rates, P)
cat( mean(hits) - 2*sd(hits)/sqrt(N), mean(hits), mean(hits) + 2*sd(hits)/sqrt(N), "\n" )
