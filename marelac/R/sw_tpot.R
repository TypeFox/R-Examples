## -----------------------------------------------------------------------------
## Potential Temperature of Seawater
## -----------------------------------------------------------------------------

sw_tpot<- function (S = 35, t = 25, p, pref = 0) {
  if (any (S<0))
    stop ("Salinity should be >= 0")

  P    <- p
  Pref <- max(0, pref)
  H  <- 10*(Pref-P)
  XK <- H*sw_adtgrad(t = t, S = S, p = P)/10

  t  <- t + 0.5*XK
  Q  <- XK
  P  <- P + 0.05*H
  XK <- H*sw_adtgrad(t = t, S = S, p = P)/10

  t  <- t + 0.29289322*(XK-Q)
  Q  <- 0.58578644*XK + 0.121320344*Q
  XK <- H*sw_adtgrad(t = t, S = S, p = P)/10

  t  <- t + 1.707106781*(XK-Q)
  Q  <- 3.414213562*XK - 4.121320344*Q
  P  <- P + 0.05*H
  XK <- H*sw_adtgrad(t = t, S = S, p = P)/10

  return (t + (XK - 2.0 * Q) / 6.0)
}
