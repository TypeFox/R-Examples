## -----------------------------------------------------------------------------
## Saturated Concentration of Oxygen in Seawater
## -----------------------------------------------------------------------------

gas_O2sat <- function(S = 35, t = 25, masl = 0,
                      method = c("Weiss", "APHA", "Paul")) {

  if (any (S < 0))
    stop ("Salinity should be >= 0")
  K <-  t + 273.15
  log10 <- function(x) log(x)/log(10)

  method <- match.arg(method)
  if (any(S != 0) & (method != "Weiss"))
      warning("Salinity value ignored by this method!")
  if (any(masl != 0) & (method != "Paul"))
      warning("Sea level height ignored by this method;\n please use gas_satconc with argument P.")
  if (! checkVecLength(list(S, t, masl)))
      warning("Arguments 'S' and 't' should have the same length or length 1.")

  ret <- switch(method,
    ## American Public Health Association
    APHA  = exp(-139.34411 + (157570.1/K) - (66423080/K^2) + (12438000000/K^3)-
           (862194900000/K^4)),
    ## Weiss, R. (1970). "The solubility of nitrogen, oxygen, and argon
    ## in water and seawater". Deep-Sea Res. 17: 721-35.
    Weiss = 1.4276 * exp(-173.4292 + 249.6339 * 100 / K +
            143.3483 * log(K / 100) - 21.8492 * K / 100 +
            S * (-0.033096 + 0.014259 * K / 100 + - 0.001700 * (K / 100)^2)),
    ## Paul, L. simple approximation that respects height above see level
    Paul = (1012-0.12 * masl)/1013 * (14.674 - 13.644 * log10(1 + t/12.8))
  )
  ret
}

