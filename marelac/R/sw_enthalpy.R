## -----------------------------------------------------------------------------
## Specific Enthalpy of Seawater
## -----------------------------------------------------------------------------

sw_enthalpy <- function (S = 35, t = 25, p = P-1.013253, P = 1.013253)  {
  if (any (S<0))
    stop ("Salinity should be >= 0")

  sw_gibbs(S, t, p, dS = 0, dt = 0, dp = 0) -
    (t + 273.15) * sw_gibbs(S, t, p, dS = 0, dt = 1, dp = 0)
}
