
# Computation of CUSUM k (reference value) for given in-control ARL and threshold h (mean monitoring)

xcusum.crit.L0h <- function(L0, h, hs=0, sided="one", r=30, L0.eps=1e-6, k.eps=1e-8) {
  h.max <- xcusum.crit(0, L0, 0)
  if ( h.max < h ) stop("h too large  or  L0 far too small")
  k1 <- 0
  L0_1 <- 0
  while ( L0_1 < L0 ) {
    k1 <- k1 + .1
    L0_1 <- xcusum.arl(k1, h, 0, hs=hs, sided=sided, r=r)
  } 
  while ( L0_1 > L0 & k1 > 0.01) {
    k1 <- k1 - .01
    L0_1 <- xcusum.arl(k1, h, 0, hs=hs, sided=sided, r=r)
  }
  k2 <- k1 + .01
  L0_2 <- xcusum.arl(k2, h, 0, hs=hs, r=r)
  dk <- 1
  while ( abs(L0-L0_2) > L0.eps  &  abs(dk) > k.eps  ) {
    k3 <- k1 + ( L0 - L0_1 ) / ( L0_2 - L0_1 ) * ( k2 - k1 )
    L0_3 <- xcusum.arl(k3, h, 0, hs=hs, sided=sided, r=r) 
    # secant rule
    dk <- k3-k2
    k1 <- k2
    L0_1 <- L0_2
    k2 <- k3
    L0_2 <- L0_3
  }
  k3
}
