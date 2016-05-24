predict_s <- function(y, b.hat, X) {
  lin.pred <- as.matrix(X) %*% b.hat
  scale <- unlink_s(y, lin.pred)
  return(scale)
}

