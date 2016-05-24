##
## runit-sbx.r - SBX crossover
##
## These checks may fail sometimes! They are simply empirical checks of the
## probability of crossover.
##

N <- 10000L
f <- sbx_operator(2, 0.8, -2, 2)
parents <- matrix(c(0, 1), ncol=2)
x <- replicate(N, c(f(parents)))[1,]

test.sbxP1 <- function() {
  p <- mean(x != 0)
  message("P1 = ", p)
  checkTrue(p > 0.78 && p < 0.82)
}

test.sbxInBounds <- function() {
  checkTrue(all(x >= -2))
  checkTrue(all(x <= 2))
}
