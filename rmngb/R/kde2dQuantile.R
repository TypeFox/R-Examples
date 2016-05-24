kde2dQuantile <- function(d, X, Y, probs = .05, ...) {
  xInd <- sapply(X, function(x) whichClosest(d$x, x))
  yInd <- sapply(Y, function(x) whichClosest(d$y, x))
  zValues <- d$z[cbind(xInd, yInd)]
  quantile(zValues, probs=probs, ...)
}
