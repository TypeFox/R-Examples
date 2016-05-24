Sjll2 <- function(b.hat, X, y, p, offset = 0, ...) {
  y.hat <- predict(y, b.hat[1:p], X[,1:p], offset)
  scale.hat <- predict_s(y, b.hat[-(1:p)], X[,-(1:p)])
  sum(jll2(y, y.hat, scale.hat, ...))
}
