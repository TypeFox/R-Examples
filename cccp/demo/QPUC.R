##
## Demo for solving an unconstrained QP
##
## Creating objects for QP
n <- 4L
M <- matrix(rnorm(n^2), nrow = n, ncol = n)
P <- crossprod(M)
q <- rnorm(n)
## Solving QP
ans <- cccp(P = P, q = q)
ans
getx(ans)
