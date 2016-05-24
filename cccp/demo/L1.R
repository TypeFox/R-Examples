##
## Demo for solving a L1-norm approximation by means of a Linear Programs
## (Example taken from cvxopt's userguide)
##
## Creating Problem
set.seed(12345)
n <- 10
m <- 20
P <- matrix(rnorm(m * n), nrow = m, ncol = n)
q <- rnorm(m)
## Solving problem by calling wrapper-function to LP
ans <- l1norm(P = P, q = q, optctrl = ctrl())
ans
getx(ans)[1:n]
