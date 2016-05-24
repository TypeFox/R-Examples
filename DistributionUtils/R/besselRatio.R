besselRatio <- function(x, nu, orderDiff, useExpScaled = 700) {
  if (x > useExpScaled) {
    besselK(x, nu + orderDiff, expon.scaled = TRUE) /
    besselK(x, nu, expon.scaled = TRUE)
  } else {
    besselK(x, nu + orderDiff) / besselK(x, nu)
  }
}
