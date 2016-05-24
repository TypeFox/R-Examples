predict.expFamily <- function(object, b.hat, X, offset = 0, m, a, ...) {
  lin.pred <- as.matrix(X) %*% b.hat + offset
  y.hat <- unlink(object, lin.pred, m, a)
  return(y.hat)
}
