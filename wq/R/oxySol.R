oxySol <-
function (t, S, P = NULL) {

  T = t + 273.15  # deg K
  lnCstar = -139.34411 + 157570.1/T - 66423080/T^2 + 1.2438e+10/T^3 -
    862194900000/T^4 - S * (0.017674 - 10.754/T + 2140.7/T^2)
  Cstar1 <- exp(lnCstar)
  if (is.null(P)) {
    # equilibrium DO Cstar at P = 1 atm
    Cstar1
  } else {   
    # transform for nonstandard pressure
    Pwv = (1 - 0.000537 * S) * exp(18.1973 * (1 - 373.16/T) +
            3.1813e-07 * (1 - exp(26.1205 * (1 - T/373.16))) - 0.018726 *
            (1 - exp(8.03945 * (1 - 373.16/T))) + 5.02802 * log(373.16/T))
    theta = 0.000975 - 1.426e-05 * t + 6.436e-08 * t^2
    Cstar1 * P * (1 - Pwv/P) * (1 - theta * P)/((1 - Pwv) * (1 - theta))
  }
}
