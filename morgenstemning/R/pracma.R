.pchip <- function(xi, yi, x) {
  h <- diff(xi)
  delta <- diff(yi)/h
  d <- .pchipslopes(h, delta)
  n <- length(xi)
  a <- (3 * delta - 2 * d[1:(n-1)] - d[2:n]) / h
  b <- (d[1:(n-1)] - 2 * delta + d[2:n])/h^2
  k <- rep(1, length(x))
  for (j in 2:(n-1)) {
    k[xi[j] <= x] <- j
  }
  s <- x - xi[k]
  v = yi[k] + s * (d[k] + s * (a[k] + s * b[k]))
  return(v)
}

.pchipslopes <- function (h, delta) {
  n <- length(h) + 1
  d <- numeric(length(h))
  k <- which(sign(delta[1:(n-2)]) * sign(delta[2:(n-1)]) > 0) + 1
  w1 <- 2 * h[k] + h[k-1]
  w2 <- h[k] + 2 * h[k-1]
  d[k] <- (w1 + w2) / (w1 / delta[k-1] + w2 / delta[k])
  d[1] <- .pchipend(h[1], h[2], delta[1], delta[2])
  d[n] <- .pchipend(h[n-1], h[n-2], delta[n-1], delta[n-2])
  return(d)
}

.pchipend <- function (h1, h2, del1, del2) {
  d <- ((2 * h1 + h2) * del1 - h1 * del2) / (h1 + h2)
  if (sign(d) != sign(del1)) {
    d <- 0
  } else if ((sign(del1) != sign(del2)) && (abs(d) > abs(3 * del1))) {
    d <- 3 * del1
  }
  return(d)
}
