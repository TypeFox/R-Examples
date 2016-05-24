#' Negative 2 log likelihood
#'
#' Calculate the -2 * log likelihood of a dataset given a specified model.
#'
#' @param b intercept and coefficients of a generalized linear model.
#' @param dataset a test dataset used to derive the likelihood.
#' @return the function returns the -2 * log likelihood.
#' @examples
#'## Using the mtcars dataset
#'## Resample, fit an ordinary least squares model and calculate likelihood
#' data(mtcars)
#' mtc.data <- cbind(1,datashape(mtcars, y = 8, x = c(1, 6, 9)))
#' head(mtc.data)
#' mtc.boot <- randboot(mtc.data, replace = TRUE)
#' boot.betas <- ml.rgr(mtc.boot)
#' loglikelihood(b = boot.betas, dataset = mtc.data)


loglikelihood <- function(b, dataset){

  dataset <- as.matrix(dataset)
  b <- matrix(b, ncol = 1)
  nc <- dim(dataset)[2]  # nc is the number of columns in d
  y <- dataset[, nc]      # y is the last column in d (a vector)
  n <- length(y)   # n is the number of rows in d (scalar)
  DS <- dataset
  DS <- DS[, 1:nc - 1]    # remove the last column from D


  linear <- DS %*% b   # linear is the matrix product of D and b

  z <- t(y) %*% linear - rep(1, n) %*% log(1 + exp(linear))
  z <- -2 * z

  return(z)

}


