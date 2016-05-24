library(PEIP)
set.seed(11)

n <- 5
A <- matrix(runif(n*n),nrow=n)
B <- matrix(runif(n*n),nrow=n)

GSVD(A, B)
