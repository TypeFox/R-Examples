MAXError <-
function(actual, prediction) {
  if (length(actual) != length(prediction)) stop("actual and prediction have different lengths")
    
  res <- max(abs(actual-prediction))
  res
}
