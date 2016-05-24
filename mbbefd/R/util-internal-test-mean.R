#internal test function of expectation for one-inflated distribution

tmean1 <- function(doifun, ...)
  integrate(function(x) x*doifun(x, ...), 0, 1)$value + 1*doifun(x=1, ...)
tmean2 <- function(poifun, ...)
  integrate(function(x) 1-poifun(x, ...), 0, 1)$value
tmean3 <- function(ecoifun, eps=1e-9, ...)
  eps/ecoifun(eps, ...)

