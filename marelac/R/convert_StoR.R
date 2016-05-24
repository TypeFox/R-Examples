## -----------------------------------------------------------------------------
## Compute Conductivity Ratio from Salinity, Temperature, and Pressure
## -----------------------------------------------------------------------------

convert_StoR <- function(S = 35, t = 25, p = max(0, P-1.013253), P = 1.013253) {
  if (any (S<0))
    stop ("Salinity should be >= 0")
  fun <- function(x) convert_RtoS(x, t = t, p = p) - S
  cond <- uniroot(fun, c(0, 5), tol = 1e-10)$root
  return(cond)
}

