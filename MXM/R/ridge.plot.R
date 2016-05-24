################################
#### Univariate ridge regression
#### Plot to see how things go
#### How the coefficients shrink
################################

#### usage: ridge.plot(target, dataset, lambda = seq(0, 5, by = 0.1) ) 

ridge.plot <- function(target, dataset, lambda = seq(0, 5, by = 0.1) ) {
  ## target is the dependent variable and can only be a vector
  ## dataset contains the independent, CONTINUOUS ONLY, variables
  ## lambda contains a grid of values of the ridge regularization parameter
  target <- as.vector(target)
  dataset <- as.matrix(dataset)
  n <- length(target)  ## sample size
  p <- ncol(dataset)  ## dimensionality of dataset
  R <- length(lambda)
  be <- matrix(nrow = p, ncol = R)
  yy <- target - mean(target)  ## center the dependent variables
  xx <- scale(dataset)[1:n, ]  ## standardize the independent variables
  sa <- svd(xx)
  u <- sa$u ;  d <- sa$d  ;  v <- sa$v
  for (i in 1:R) {
    be[, i] <- ( v %*% diag( d/( d^2 + lambda[i] ) ) %*% t(u) ) %*% yy 
  }
  plot(lambda, be[1,], type = "l", col = 1, lty = 1,
  ylim = c( min(be), max(be) ), xlab = expression(paste(lambda, " values")), 
  ylab = "Beta coefficients")
  for (i in 2:p)  lines(lambda, be[i, ], col = i, lty = i) 
}