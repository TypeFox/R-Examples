missingentirerow <- function(x) {
  n <- length(x)
  nm <- sum(is.na(x))
  y <- n==nm
  return(y)
}
