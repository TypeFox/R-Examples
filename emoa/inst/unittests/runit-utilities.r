##
## runit-utilities.r - check internal utility functions
##

test.inbounds <- function() {
  x <- c(0, 1, 2)
  checkEquals(emoa:::inbounds(x, 0, 1), c(0, 1, 1))
  checkEquals(emoa:::inbounds(x, 1, 2), c(1, 1, 2))
  checkEquals(emoa:::inbounds(x, 1, 1), c(1, 1, 1))

  checkEquals(emoa:::inbounds(x, c(0, 0, 1), c(1, 1, 2)), c(0, 1, 2))
  checkEquals(emoa:::inbounds(x, c(1, 0, 1), c(1, 1, 2)), c(1, 1, 2))
  checkEquals(emoa:::inbounds(x, c(0, 0, 1), c(1, 1, 1)), c(0, 1, 1))
}
