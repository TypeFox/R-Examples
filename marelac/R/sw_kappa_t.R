## -----------------------------------------------------------------------------
## Isothermal Compressibility of Seawater
## -----------------------------------------------------------------------------

sw_kappa_t <- function (S = 35, t = 25, p = P-1.013253, P = 1.013253)  {
  if (any (S < 0))
    stop ("Salinity should be >= 0")

  -1e4 * sw_gibbs(S, t, p, dS = 0, dt = 0, dp = 2) /
      sw_gibbs(S, t, p, dS = 0, dt = 0, dp = 1)
}
