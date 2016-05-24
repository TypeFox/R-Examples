cubic.spline.penalty2 <-
function(knots, interval)
{
  tau <- c(interval[1L], knots, interval[2L])
  k <- length(tau)
  i1 <- 1L : (k - 1L)
  i2 <- 2L : k
  tau.mid <- (tau[i1] + tau[i2]) / 2
  Z <- cubic.spline.deriv2(c(tau, tau.mid), knots, interval)
  i3 <- (k + 1L) : (2L * k - 1L)
  dtau <- tau[i2] - tau[i1]
  
  # uses Simpson's rule for numerical integration, whose result
  # is actually exact since the integrand is at most quadratic
  D1 <- crossprod(Z[i1, ], Z[i1, ] * dtau) / 6
  D2 <- crossprod(Z[i2, ], Z[i2, ] * dtau) / 6
  D3 <- crossprod(Z[i3, ], Z[i3, ] * dtau) * 4 / 6
  D1 + D2 + D3
}
