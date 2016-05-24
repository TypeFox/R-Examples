p.unconstrain <- function(p, lower = -Inf, upper = Inf, f = 1) {
  if (!(min(lower) == -Inf | max(upper) == Inf)) {
    p <- tan(pi/2 * (2 * p - upper - lower) / (upper - lower)) / f
  }
  p
}

p.constrain <- function(p, lower = -Inf, upper = Inf, f = 1) {
  if (!(min(lower) == -Inf | max(upper) == Inf)) {
    p <- 1/2 * (upper + lower) + (upper - lower) * atan(f * p)/pi
  }
  p
}
