`aweibull6` <-
function(p, lower = 0, upper = 365) {
  integrate(f = fweibull6, lower = lower, upper = upper, p = p)$value
}

