`print.anacor` <-
function(x, ...)
{
# x ... object of class "anacor"
  cat("\nCA fit: \n")
  cat("Sum of eigenvalues: ", sum(x$eigen.values),"\n")
  if (!is.null(x$rmse[[1]])) cat("Benzecri RMSE rows: ", x$rmse[[1]],"\n")
  if (!is.null(x$rmse[[2]])) cat("Benzecri RMSE columns: ", x$rmse[[2]],"\n")
  cat("\nTotal chi-square value:",round(sum(x$chisq.decomp[,1]),3),"\n")
  cat("\nChi-Square decomposition: \n")
  print(round(x$chisq.decomp,3))
  cat("\n")
}

