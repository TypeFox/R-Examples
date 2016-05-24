dgp <- function(n, b, intercept = TRUE, sc = -1) {
  X <- as.matrix(rep(1, times = n))
  if (intercept == FALSE) {
  X <- X[, -1]
  }
for (i in 1:(length(b)-1)) {
  x <- runif(n, 0, 50)
  X <- cbind (X, x)
  }
  x <- X
  if (intercept) {
      x <- X[, -1]
      }
y_ols <- X%*%b
u <- abs(rnorm(n, 0, 1))
#u <- rhalfnorm(n, 1)
v <- rnorm(n, 0, 1)
epsilon <- v - sc*u
y <- y_ols + epsilon
return(list(X = X, x = x, y = y, u = u, v = v))
}