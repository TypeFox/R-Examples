set.seed(779)
`%perp%` <- function(y,x)
      lsfit(x,y, intercept  = FALSE)$residuals
x <- 1:10
y <- x + .25 * rnorm(10)
y %perp% x
