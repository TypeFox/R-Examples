## -----------------------------------------------------------------------------
## Isentropic Compressibility of Seawater
## -----------------------------------------------------------------------------

sw_kappa <- function (S = 35, t = 25, p = P-1.013253, P = 1.013253) {

  if (any (S < 0))
    stop ("Salinity should be >= 0")

  g_tt <- sw_gibbs(S, t, p, dS = 0, dt = 2, dp = 0)
  g_tp <- sw_gibbs(S, t, p, dS = 0, dt = 1, dp = 1)

  1e4 * (g_tp * g_tp - g_tt * sw_gibbs(S, t, p, dS = 0, dt = 0, dp = 2)) /
    (sw_gibbs(S, t, p, dS = 0, dt = 0, dp = 1) * g_tt)
}
