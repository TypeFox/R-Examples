bspline.deriv <-
function(x, order, knots, interval, q = 0L)
{
  q <- as.integer(q)
  if (q <= 0L) {
    return(bspline(x, order, knots, interval))
  }

  p <- length(knots) # number of interior knots
  m <- as.integer(order) # m == 4 means cubic spline
  tau <- c(rep(interval[1], m), knots, rep(interval[2], m))
  
  Z <- rbind(0, t(bspline.deriv(x, m - 1L, knots, interval, q - 1L)), 0)
  k <- m - 1L # degree of the polynomial 
  i <- 1L : (m + p)
  a1 <- tau[i + k] - tau[i]
  a1 <- ifelse(a1 < 1e-12, 0, 1 / a1)
  a2 <- tau[i + k + 1L] - tau[i + 1L]
  a2 <- ifelse(a2 < 1e-12, 0, 1 / a2)
  Z <- (k * a1) * Z[i, ] - (k * a2) * Z[i + 1L, ]
  
  t(Z)
}
