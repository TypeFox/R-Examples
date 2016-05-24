bspline <-
function(x, order, knots, interval) 
{
  p <- length(knots) # number of interior knots
  m <- as.integer(order) # m == 4 means cubic spline
  tau <- c(rep(interval[1L], m), knots, rep(interval[2L], m))
  
  Z <- matrix(0, 2L * m + p, length(x))
  i <- m : (m + p - 1L)
  Z[i, ] <- outer(tau[i], x, "<=") * outer(tau[i + 1L], x, ">")
  i <- m + p
  Z[i, ] <- outer(tau[i], x, "<=") * outer(tau[i + 1L], x, ">=")
  
  k <- 1L
  while (k < m) {
    i <- (m - k) : (m + p)
    a1 <- tau[i + k] - tau[i]
    a1 <- ifelse(a1 < 1e-12, 0, 1 / a1)
    a2 <- tau[i + k + 1L] - tau[i + 1L]
    a2 <- ifelse(a2 < 1e-12, 0, 1 / a2)
    Z[i, ] <- -a1 * outer(tau[i], x, "-") * Z[i, ] + 
      a2 * outer(tau[i + k + 1L], x, "-") * Z[i + 1L, ]
    k <- k + 1L
  }
  
  t(Z[1L : (m + p), ])
}
