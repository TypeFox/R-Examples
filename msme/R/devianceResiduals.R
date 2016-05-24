devianceResiduals <- function(y, b.hat, X, offset = 0, ...) {
  y.hat <- predict(y, b.hat, X, offset)
  sign(y - y.hat) * sqrt(2 * (jll(y, y, ...) - jll(y, y.hat, ...)))
}
