## This script is purely experimental
## and contains heuristics for initializing start values

init_twostep <- function(x, y) {
  ## simple heuristics
  m <- lm(y~x)
  yi <- ya <- unname(min(y) / 2)
  K  <- unname(max(y))
  mumax <- unname(coef(m)[2])
  kw <- unname(0.1)

  c(yi=yi, ya=ya, kw=kw, mumax=mumax, K=K)
}

init_logistic <- function(x, y) {
  ## simple heuristics
  m  <- lm(y~x)
  y0 <- unname(min(y))
  K  <- unname(max(y))
  mumax <- unname(coef(m)[2])

  c(y0=y0, mumax=mumax, K=K)
}

## data(bactgrowth)
## p <- init_twostep(dat$time, dat$value)
