library(scs)
m <- 100
n <- 100

A <- matrix(rnorm(m*n),m,n)
b <- rnorm(m)
objective = rnorm(n)
cone = list(f = m)
params = list(eps = 1e-3, max_iters = 5000)
sol <- scs(A, b, objective, cone, params)
sol

