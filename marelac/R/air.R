## -----------------------------------------------------------------------------
## Air density and specific humidity
## -----------------------------------------------------------------------------

air_spechum <- function (t = 25, rh = 50, P = 1.013253) {
  if (! checkVecLength(list(t, rh, P)))
    warning("Arguments 't', 'rh' and 'P' should have the same length or length 1.")

  MH20 <- 18.01534
  Mdry <- 28.9644
  p    <- P*1e3
  e    <- vapor.hPa(t = t)
  PH2O <- rh*e/100
  xH2O <- PH2O/p
# specific humidity

  xH2O*MH20/(xH2O*MH20+(1-xH2O)*Mdry)
}


air_density <- function (t = 25, P = 1.013253) {
  if (! checkVecLength(list(t, P)))
    warning("Arguments 't' and 'P' should have the same length or length 1.")
  Mdry <- 28.9644  # g/mol
  p <- P*1e5
  R <- 8.31447215
  p / (R *(t + 273.15)) * Mdry / 1000
}
