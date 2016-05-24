rmse <- function(forecast, validation) {
  stopifnot(is.vector(forecast, mode='numeric'))
  stopifnot(is.vector(validation, mode='numeric'))
  stopifnot(!anyNA(forecast))
  stopifnot(!anyNA(validation))
  stopifnot(length(forecast) == length(validation))
  #stopifnot(length(history) > 1)

  return(sqrt(mean((forecast - validation)^2)))
}
