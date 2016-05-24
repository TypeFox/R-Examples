mase <- function(forecast, validation) {
  stopifnot(is.vector(forecast, mode='numeric'))
  stopifnot(is.vector(validation, mode='numeric'))
  stopifnot(!anyNA(forecast))
  stopifnot(!anyNA(validation))
  stopifnot(length(forecast) == length(validation))
  #stopifnot(length(history) > 1)

  n <- length(validation)

  #mase <- mean(abs(validation - forecast)) / mean(abs(history[-length(history)] - history[-1]))
  mase <- sum(abs(validation - forecast)) / sum(abs(validation[-1] - validation[-n])) / (n/(n-1))

  return(mase)
}
