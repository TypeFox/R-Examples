lower_triangle <- function (x) {
  n <- nrow(x)
  return(x[outer (1:n, 1:n, ">" )])
}

fill_symmetric <- function (x) {
  m <- length (x)
  n <- 0.5 + sqrt (0.25 + 2 * m)
  d <- matrix (0, n, n)
  d[outer (1:n, 1:n, ">" )] <- x
  return(d + t(d))
}

Euclid <- function (x) {
  c <- tcrossprod (x)
  d <- diag (c)
  return(sqrt (outer (d, d, "+" ) - 2 * c))
}

circular <- function (n) {
  x <- seq (0, 2 * pi, length = n + 1)
  z <- matrix (0, n + 1, 2)  
  z[, 1] <- sin (x)
  z[, 2] <- cos (x)
  return (z[-1, ])
}

direct_sum <- function (x) {
  n <- length (x)
  nr <- sapply (x, nrow)
  nc <- sapply (x, ncol)
  s <- matrix (0, sum (nr), sum (nc))
  k <- 0
  l <- 0
  for (j in 1 : n) {
    s[k + (1 : nr[j]), l + (1 : nc[j])] <- x[[j]]
    k <- k + nr[j]
    l <- l + nc[j]
  }
  return (s)
}
