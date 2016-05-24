#' Sum of Square Errors
#'
#' Compute the sum of squared prediction errors (or residual sum of squares)
#' when a linear model is applied to a dataset.
#'
#' @param b vector or column-matrix of regression coefficients
#' @param dataset a matrix or dataframe. The final column is the outcome variable.
#' @return The function returns the sum of square errors.
#'
#' @examples
#'## Using simulated data derived from the iris dataset
#' mu <- c(rep(0, 4))
#' covmatr <- matrix(c(0.7, -0.04, 1.3, 0.5, -0.04, 0.2, -0.3, -0.1,
#' 1.3, -0.3, 3.1, 1.3, 0.5, -0.1, 1.3, 0.6), ncol = 4)
#' sim.dat <- randnor(n = 100, mu = mu, Cov = covmatr)
#' sim.dat <- cbind(1, sim.dat)

#'## resample and fit an ordinary least squares model, and then
#'## calculate the sum of square errors of the model when applied
#'## to the original data
#' sim.boot <- randboot(sim.dat, replace = TRUE)
#' boot.betas <- ols.rgr(sim.boot)
#' sse(b = boot.betas, dataset = sim.dat)

sse <- function(b, dataset) {

  b <- matrix(b, ncol = 1)
  dataset <- as.matrix(dataset)
  if (dim(b)[1] != (dim(dataset)[2] - 1)) stop("The number of regression
                   coefficients does not match the number of predictor variables")

  rd <- dim(dataset)[1]
  d0 <- dataset[, (1:dim(dataset)[2] - 1)]
  yhat <- d0 %*% b
  res <- dataset[, dim(dataset)[2]] - yhat
  z <- t(res) %*% res
  return(z)
}
