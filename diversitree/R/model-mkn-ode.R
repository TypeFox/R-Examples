## ODE interface to MK models.

initial.tip.mkn.ode <- function(cache) {
  k <- cache$info$k
  y <- matrix(0, k+1, k, TRUE)
  y[k+1,] <- diag(y[1:k,]) <- 1
  y <- matrix.to.list(y)
  y.i <- cache$states
  dt.tips.grouped(y, y.i, cache)
}

make.branches.mkn.ode <- function(cache, control)
  make.branches.dtlik(cache$info, control)

derivs.mkn.ode <- function(t, y, pars) {
  k <- length(y)
  Q <- matrix(pars, k, k)
  ret <- Q %*% y
  dim(ret) <- NULL
  ret
}
