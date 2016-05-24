MAPE <-
function(actual, prediction) {
  if (length(actual) != length(prediction)) stop("actual and prediction have different lengths")
  
  n <- length(actual)
  
  res <- (100 / n) * sum(abs((actual-prediction)/actual))
  res
}
