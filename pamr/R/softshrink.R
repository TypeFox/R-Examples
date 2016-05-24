soft.shrink <-function(delta, threshold) {
  dif <- abs(delta) - threshold
  delta <- sign(delta) * dif * (dif > 0)
  nonzero <- sum(drop((dif > 0) %*% rep(1, ncol(delta))) > 0)
  attr(delta, "nonzero") <- nonzero
  delta
}
