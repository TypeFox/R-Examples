z.score <- function (k, n, p=.5, correct=TRUE) {
  if (any(k < 0) || any(k > n) || any(n < 1)) stop("arguments must be integer vectors with 0 <= k <= n")
  if (any(p < 0) || any(p > 1)) stop("null hypothesis proportion p must be in range [0,1]")

  l <- max(length(k), length(n), length(p)) # ensure that all vectors have the same length
  if (length(k) < l) k <- rep(k, length.out=l)
  if (length(n) < l) n <- rep(n, length.out=l)
  if (length(p) < l) p <- rep(p, length.out=l)

  expected <- n * p                     # compute z-score (with optional Yates' correction)
  variance <- n * p * (1-p)
  d <- k - expected
  if (correct) d <- sign(d) * pmax(0, abs(d) - .5)
  z <- d / sqrt(variance)
  z
}
