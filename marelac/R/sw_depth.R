## -----------------------------------------------------------------------------
## Water Depth from Hydrostatic Pressure and Latitude
## -----------------------------------------------------------------------------

sw_depth <- function (p = P-1.013253, P = 1.013253, lat = 0) {
  P <- p*10    # P=hydrostatic pressure, in dbar
  denom <- gravity(lat) + 1.092e-6*P
  nom   <- (9.72659 + (-2.2512e-5 + (2.279e-10 - 1.82e-15*P)*P)*P)*P

  return (nom / denom)
}
