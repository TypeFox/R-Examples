cubic.spline.penalty1 <-
function(knots, interval)
{
  tau <- c(interval[1L], knots, interval[2L])
  k <- length(tau)
  i1 <- 1L : (k - 1L)
  i2 <- 2L : k
  tau.mid0 <- tau[i1]
  tau.mid4 <- tau[i2]
  tau.mid2 <- (tau.mid0 + tau.mid4) / 2
  tau.mid1 <- (tau.mid0 + tau.mid2) / 2
  tau.mid3 <- (tau.mid2 + tau.mid4) / 2
  Z <- cubic.spline.deriv1(c(tau.mid0, tau.mid1, tau.mid2, 
    tau.mid3, tau.mid4), knots, interval)
  i <- list("0" = i1, "1" = i1 + (k - 1L), "2" = i1 + (k - 1L) * 2L,
    "3" = i1 + (k - 1L) * 3L, "4" = i1 + (k - 1L) * 4L)
  dtau <- tau[i2] - tau[i1]
  
  # uses Boole's rule for numerical integration, whose result
  # is actually exact since the integrand is at most of degree 4
  D0 <- crossprod(Z[i[["0"]], ], Z[i[["0"]], ] * dtau) * 7 / 90
  D1 <- crossprod(Z[i[["1"]], ], Z[i[["1"]], ] * dtau) * 32 / 90
  D2 <- crossprod(Z[i[["2"]], ], Z[i[["2"]], ] * dtau) * 12 / 90
  D3 <- crossprod(Z[i[["3"]], ], Z[i[["3"]], ] * dtau) * 32 / 90
  D4 <- crossprod(Z[i[["4"]], ], Z[i[["4"]], ] * dtau) * 7 / 90
  D0 + D1 + D2 + D3 + D4
}
