library(scs)

A <- matrix(c(1, 1), ncol=1)
b <- c(1, 1)
c <- c(1)
cone <- list(f = 2)
params <- list(eps = 1e-3, max_iters = 50, verbose=TRUE)
sol <- scs(A, b, c, cone, params) 
sol


