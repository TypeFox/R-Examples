## -----------------------------------------------------------------------------
## Conversion Between Different Barometric Units
## -----------------------------------------------------------------------------

convert_p <- function(x, unit = c("Pa", "bar", "at", "atm", "torr")) {
  if (!is.numeric(x)) stop("x must be numeric")
  ## factors taken from Wikipedia
  unit    <- match.arg(unit)
  factors <- c(Pa = 1, bar = 1e-5, at = 1.0197e-5, atm = 9.8692e-6, torr = 7.5006e-3)
  units   <- names(factors)
  p       <- match(unit, units)
  as.data.frame(x %o% factors / factors[p])
}

