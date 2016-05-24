smape <- function(forecast, validation) {
  stopifnot(is.vector(forecast, mode='numeric'))
  stopifnot(is.vector(validation, mode='numeric'))
  stopifnot(!anyNA(forecast))
  stopifnot(!anyNA(validation))
  stopifnot(length(forecast) == length(validation))

  r <- abs(forecast - validation) / ((abs(forecast) + abs(validation)) / 2)
  r <- mean(r)
  return(r)
}
