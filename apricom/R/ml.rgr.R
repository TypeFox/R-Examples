#' Logistic Regression with Maximum Likelihood Estimation
#'
#' Fit a logistic regression model using maximum likelihood estimation
#'
#' This function is a wrapper for \code{\link[stats]{glm.fit}}, for convenient
#' application within several functions in the \pkg{apricomp} package. This function may
#' be called directly. For regression with an intercept included, the first column in
#' the dataset must be a column of 1s.
#'
#' @importFrom stats glm.fit binomial
#'
#' @param dataset a \code{p} x \code{m} data matrix, where the final column is a
#'        binary outcome variable. \code{datashape} may be applied to data so that
#'        the dataset is in the correct format for this function (see manual)
#' @return The function returns a column-vector containing the logistic
#'         regression coefficients and intercept (if specified).
#'
#' @examples
#'## Logistic regression using a subset of the mtcars data (outcome is "vs")
#' data(mtcars)
#' mtc.df <- mtcars[, c(8, 1, 9)]
#' mtc.shaped <- datashape(dataset = mtc.df, y = 1)
#' ml.rgr(mtc.shaped)
#' ml.rgr(cbind(1,mtc.shaped))

ml.rgr <- function(dataset) {

  nc <- dim(dataset)[2]
  X <- dataset[, 1:nc - 1]
  Y <- dataset[, nc]

# fit logistic regression model using glm.fit
  ml.fit <- glm.fit(X, Y, family = binomial(link = "logit"))
  coeffs <- as.matrix(ml.fit$coefficients, ncol = 1)

  return(coeffs)

}
