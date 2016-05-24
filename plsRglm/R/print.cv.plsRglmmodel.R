print.cv.plsRglmmodel <- function(x,...)
{
  cat("Number of repeated crossvalidations:\n")
  print(length(x$results_kfolds))
  cat("Number of folds for each crossvalidation:\n")
  print(length(x$results_kfolds[[1]]))
}
