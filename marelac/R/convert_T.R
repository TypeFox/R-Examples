## -----------------------------------------------------------------------------
## Conversion Between Different Temperature Units
## -----------------------------------------------------------------------------

convert_T <- function(x, unit = c("K", "C", "F")) {
  if (!is.numeric(x)) stop("x must be numeric")
  unit   <- match.arg(unit)
  const  <- c(K = 0, C = -273.15, F = -459.67)
  a      <- c(K = 1, C = 1,       F = 5/9)
  units <- names(const)
  u <- match(unit, units)
  ## 1) convert everything to K
  TK <- (x - const[u]) * a[u]
  ## 2) convert to the other scales
  ret <- data.frame(K = TK, C = TK - 273.15, F = TK * 1.8 - 459.67)
  row.names(ret) <- NULL
  ret
}

