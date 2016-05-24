var.asymp <- function(p, n = 1) {
  p * (1 - p)/n
}

var.cloglog <- function(p, n = 1) {
  mu <- log(p)
  (1 - p)/n/p/mu^2
}

var.logit <- function(p, n = 1) {
  1/n/p/(1 - p)
}

var.probit <- function(p, n = 1) {
  z <- qnorm(p)
  p * (1 - p)/n/dnorm(z)^2
}

binom.lchoose <- function(n, x) {
  lgamma(n + 1) - lgamma(n - x + 1) - lgamma(x + 1)
}

ldbinom <- function(x, size, prob, log = TRUE) {
  log.f <- binom.lchoose(size, x) + x * log(prob) + (size - x) * log(1 - prob)
  if(log) log.f else exp(log.f)
}
