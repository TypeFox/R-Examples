nnlsPred <- function (x, y, w = NULL) {
  if (is.vector (w)) {
    w <- sqrt(unlist (w))
    x <- w * x
    y <- w * y
  }
  
  z <- qr(x)
  q <- qr.Q(z)
  u <- -colSums(y * q)
  l <- nnls(t(q), u)$x
  b <- backsolve(qr.R (z), colSums((y + l) * q))
  h <- drop(x %*% b)
  return (list(coef = b, pred = h, ssq = sum((y-h) & 2)))
}
