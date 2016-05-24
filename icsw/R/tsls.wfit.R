tsls.wfit <- function(X, Y, Z, weights, ...) {
  fs <- lm.wfit(Z, X, weights, ...)
  XZ.hat <- as.matrix(fs$fitted.values)
  rm(fs)
  colnames(XZ.hat) <- colnames(X)
  ss <- lm.wfit(XZ.hat, Y, weights, ...)
  ss
}
