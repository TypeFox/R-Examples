##
## runit-ei.r - Epsilon Indicator tests
##

x <- matrix(c(1.0, 0.5, 0.0,
              0.0, 0.5, 1.0),            
            ncol=3, byrow=TRUE)

runit.epsilon_indicator <- function() {
  k <- nrow(x)
  n <- ncol(x)
  ## Check for different permutations of the rows and columns of
  ## points.
  for (i in 1:10) {
    o <- sample(1:n)
    p <- sample(1:k)
    m <- x[p,o]
    for (delta in seq(0, 1, by=0.2)) {
      checkEqualsNumeric(epsilon_indicator(m, m + delta), -delta)
      checkEqualsNumeric(epsilon_indicator(m + delta, m), delta)
    }
  }
  ## Check different sized matrices:
  checkEqualsNumeric(epsilon_indicator(x, x[,-2] + 0.2), -0.2)
  ## Negative values:
  checkException(epsilon_indicator(x, x-10))
  checkException(epsilon_indicator(x-10, x))
}

ruint.epsilon_indicator_differently_sized_fronts <- function() {
  front1 <- matrix(runif(100), nrow=2)
  front2 <- matrix(runif(150), nrow=3)
  checkException(checkepsilon_indicator(nondominated_points(front1),
                                        nondominated_points(front2)))
}
