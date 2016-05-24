##
## runit-poly.r - Polynomial mutation
##
## These checks may fail sometimes! They are simply empirical checks of the
## probability of crossover.
##

N <- 10000L
f <- pm_operator(5, 0.8, -10, 10)
x <- replicate(N, f(5))

test.polyP1 <- function() {
  p <- mean(x != 5)
  message("P1 = ", p)
  checkTrue(p > 0.78 && p < 0.82)
}

test.polyP2 <- function() {
  p <- mean(x < 5)
  message("P2 = ", p)
  checkTrue(p > 0.38 && p < 0.42)
}

test.polyInBounds <- function() {
  checkTrue(all(x >= -10))
  checkTrue(all(x <= 10))
}
