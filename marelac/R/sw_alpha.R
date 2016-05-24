## -----------------------------------------------------------------------------
## Thermal Expansion Coefficient of Seawater
## -----------------------------------------------------------------------------

sw_alpha <- function (S = 35, t = 25, p = P-1.013253, P = 1.013253) {
  if (any (S<0))
    stop ("Salinity should be >= 0")

  sw_gibbs(S, t, p, dS = 0, dt = 1, dp = 1) /
    sw_gibbs(S, t, p, dS = 0, dt = 0, dp = 1)
}
