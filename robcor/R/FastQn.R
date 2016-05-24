.ONE.OVER.SQRT2 <- 1 / sqrt(2)

FastQn <- function(x, center = median(x), scale = mad(x, center)) {
  n <- length(x)
  u2 <- ((x - center) / scale) ^ 2
  z <- colSums(outer(u2, c(0, 1), "^") * exp(-0.5 * u2))
  scale * (1 - (z[1] - n * .ONE.OVER.SQRT2) / z[2])
}

fqn <- FastQn

s_FastQn <- function(x, mu.too = FALSE, center = median(x), ...) {
  c(if (mu.too) center, FastQn(x, center = center, ...))
}
