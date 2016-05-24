fitDistances <- function(x, family) {
  getDistances <- function(x) {
    rep(x$distances, x$seed.counts)
  }
  getEstimates <- function(distance) {
    fitdistr(distance, family)$estimate
  }
  distances <- lapply(x$transects, getDistances)
  parameter.list <- lapply(distances, getEstimates)
  parameters <- colMeans(do.call(rbind, parameter.list))
  return(parameters)
}
