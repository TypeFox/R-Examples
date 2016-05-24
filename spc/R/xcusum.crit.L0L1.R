
# Computation of CUSUM k (reference value) and threshold h for given in-control ARL L0 and out-of-control ARL L1 (mean monitoring)
# Ewan & Kemp 1960 or Kemp 1962

xcusum.crit.L0L1 <- function(L0, L1, hs = 0, sided="one", r = 30, L1.eps=1e-6, k.eps=1e-8) {
  k1 <- 0
  L1_1 <- L1 + 1
  while ( L1_1 > L1 ) {
    k1 <- k1 + .1
    h1   <- xcusum.crit(k1, L0, hs=hs, sided=sided, r=r)
    L1_1 <- xcusum.arl(k1, h1, 2*k1, hs=hs, sided=sided, r=r)
  } 
  while ( L1_1 < L1 & k1 > 0.01 ) {
    k1 <- k1 - .01
    h1   <- xcusum.crit(k1, L0, hs=hs, sided=sided, r=r)
    L1_1 <- xcusum.arl(k1, h1, 2*k1, hs=hs, sided=sided, r=r)
  }
  k2 <- k1 + .01
  h2   <- xcusum.crit(k2, L0, hs=hs, sided=sided, r=r)
  L1_2 <- xcusum.arl(k2, h2, 2*k2, hs=hs, sided=sided, r=r)
  dk <- 1
  while ( abs(L1-L1_2) > L1.eps  &  abs(dk) > k.eps  ) {
    k3 <- k1 + ( L1 - L1_1 ) / ( L1_2 - L1_1 ) * ( k2 - k1 )
    h3   <- xcusum.crit(k3, L0, hs=hs, sided=sided, r=r)
    L1_3 <- xcusum.arl(k3, h3, 2*k3, hs=hs, sided=sided, r=r)
    # secant rule
    dk <- k3-k2
    k1 <- k2
    L1_1 <- L1_2
    k2 <- k3
    L1_2 <- L1_3
  }
  result <- c(k3, h3)
  names(result) <- c("k", "h")
  result
}
