## -----------------------------------------------------------------------------
## Gravity from Latitude
## -----------------------------------------------------------------------------

gravity <- function (lat = 0, method = c("Moritz", "UNESCO")) {

  method <- match.arg(method)
  X <- sin(lat * pi / 180.)
  X <- X * X
  X2 <- (sin(2 * lat * pi / 180.)) ^2
  if (method == "UNESCO")
    grav <-  9.780318 * (1.0 + (5.2788e-3 + 2.36e-5 * X) *X)
  else grav <-  9.780327*(1.0 + 0.0053024 * X -0.0000058 * X2)
  return(grav)
}
