#' Checks the arguments passed to the error rate estimator functions.
#'
#' This function is a helper function that checks the arguments passed
#' to make sure the input is valid and consistent across all
#' error rate estimators.
#'
#' We expect that the first two arguments of the classifier function given in
#' \code{train} are \code{x} and \code{y}, corresponding to the data matrix and
#' the vector of their labels. Additional arguments can be passed to the
#' \code{train} function. The returned object should be a classifier that will
#' be passed to the function given in the \code{classify} argument.
#'
#' @param x a matrix of n observations and p features
#' @param y a vector of n class labels
#' @param train a function that builds the classifier (See details)
#' @param classify a function that classified observations from the constructed
#' classifier from \code{train}. (See details.)
#' @return \code{TRUE} invisibly if no errors are encountered.
check_arguments <- function(x, y, train, classify) {
  if (nrow(x) != length(y)) {
    stop("The number of observations must match the number of class labels.")
  }

  if (nlevels(y) < 2) {
    stop("There must be at least two classes given in the class labels vector 'y'.")
  }
  
  if (any(table(y) < 2)) {
    stop("Each class should have at least two observations.")
  }

  train <- try(match.fun(train), silent = TRUE)
  if (inherits(train, "try-error")) {
    stop("The 'train' function does not exist.")
  }

  if (missing(classify) || is.null(classify)) {
    classify <- match.fun("classify")
  } else {
    classify <- try(match.fun(classify), silent = TRUE)
    if (inherits(classify, "try-error")) {
      stop("The 'classify' function does not exist.")
    }
  }
  invisible(TRUE)
}
