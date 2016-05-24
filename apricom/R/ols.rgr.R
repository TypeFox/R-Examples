#' Linear Regression using Ordinary Least Squares
#'
#' Fit a linear regression model using Ordinary Least Squares.
#'
#' This function may be called directly. For regression with an intercept included,
#' the first column in the dataset must be a column of 1s.
#'
#' @param dataset a \code{p} x \code{m} data matrix, where the final column
#'         is a continuous outcome variable. \code{datashape} may be applied to data
#'         so that the dataset is in the correct format for this function (see manual)
#' @return the function returns a column-vector containing the
#'         linear regression coefficients.
#'
#' @examples
#'## Linear regression using a subset of the mtcars data (outcome is "wt")
#' data(mtcars)
#' mtc.df <- mtcars[, c(6, 1, 4)]
#' mtc.shaped <- datashape(dataset = mtc.df, y = 1)
#' ols.rgr(mtc.shaped)
#' ols.rgr(cbind(1,mtc.shaped))

ols.rgr<- function(dataset){

  dataset <- as.matrix(dataset)
  M <- dataset[, 1:(dim(dataset)[2] - 1)]
  V <- solve(t(M) %*% M)
  coeffs <- V %*% t(M) %*% dataset[, dim(dataset)[2]]
  return(coeffs)

  }
