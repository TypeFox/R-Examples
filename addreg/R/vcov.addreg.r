vcov.addreg <- function(object, ...) {
  summary.addreg(object)$cov.scaled
}