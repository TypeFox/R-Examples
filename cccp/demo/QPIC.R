##
## Demo for solving a QP with equality and inequality constraints
## (Example taken from cvxopt's userguide)
##
## Creating QP
P <- 2 * matrix(c(2, .5, .5, 1), nrow = 2, ncol = 2)
q <- c(1.0, 1.0)
G <- -diag(2)
h <- rep(0, 2)
nno1 <- nnoc(G = G, h = h)
A <- matrix(c(1.0, 1.0), nrow = 1, ncol = 2)
b <- 1.0
## Solving QP
ans <- cccp(P = P, q = q, A = A, b = b, cList = list(nno1))
ans
getx(ans)
