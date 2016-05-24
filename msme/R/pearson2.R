pearson2 <- function(object, ...) {
  dispersion <-
    with(object, {
      y.hat <- predict(y, coefficients[1:p], X[,1:p], offset)
      a.hat <- predict_s(y, coefficients[-(1:p)], X[,-(1:p)])
      pearson <- sum(((y - y.hat)^2) /
                     variance(y, y.hat, a.hat))
      pearson / df.residual})
  return(dispersion)
}
