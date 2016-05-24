sMAPE <-
function(actual, prediction) {
  if (length(actual) != length(prediction)) stop("actual and prediction have different lengths")
  
  n <- length(actual)
  
  res <- (1/n) * sum(abs(actual-prediction) / ((abs(actual)+abs(prediction))/2))
  res
}
