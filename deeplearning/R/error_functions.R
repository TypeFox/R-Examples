#' Calculates the cross entropy error
#'
#' This function calculates the cross entropy error and its first order derivatives
#'
#' @param output the output value
#' @param target the target value
#'
#' @export

crossEntropyErr <- function(output, target) {
  # err <- - sum(target[] * log(output[]) + (1 - target[]) * log(1 - output[]))
  err <- - sum(target * log(output) + (1 - target) * log(1 - output))
  err2 <- (1-target)/(1-output) - target/output

  ret <- list()
  ret[[1]] <- err
  ret[[2]] <- err2
  ret[[3]] <- "Cross Entropy Error"
  return(ret)
}

#' Calculates the mean squared error
#'
#' This function calculates the mean squared error and its first order derivatives
#'
#' @param output the output value
#' @param target the target value
#'
#' @export

meanSquareErr <- function(output, target) {
  err <- 1/2 * sum(output - target)^2 / dim(output)[[1]]
  err2 <-  (output - target)
  ret <- list()
  ret[[1]] <- err
  ret[[2]] <- err2
  ret[[3]] <- "Mean Squared Error"
  return(ret)

}

#' Calculates the classification error
#'
#' This function calculates the classification error
#'
#' @param output the output of a classifier in the form of probability. Probability > 1
#' will be treated as positive (target = 1).
#' @param target the target variable
#'
#' @export

classification_error <- function(output, target) {
  boolOut <- (output > 0.5) * 1
  boolOutTarget <- cbind(boolOut, target)
  rows <- nrow(target)
  cols <- ncol(target)
  classification_error <- sum(apply(boolOutTarget, 1, function(y)
  { any(y[1:cols] != y[(cols+1):(2*cols)])})) / rows * 100

  ret <- list()
  ret[[1]] <- classification_error
  ret[[2]] <- "Classification Error"
  return (ret)
}
