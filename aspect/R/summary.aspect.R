summary.aspect <- function(object, ...)
{
  cat("\nCorrelation matrix of the scaled data:\n")
  print(object$cormat)
  cat("\n")
  if (!is.null(object$cor.rat)) {
    cat("\nCorrelation ratios:\n")
    print(object$cor.rat)
  }
  if (!is.null(object$eigencor)) {
    cat("\nEigenvalues of the correlation matrix:\n")
    print(object$eigencor)
  }
  cat("\n")
  cat("Category scores:\n")
  namevec <- names(object$catscores)
  for (i in 1:length(object$catscores)) {
    cat(namevec[i],":\n ", sep = "")
    print(object$catscores[[i]])
    cat("\n")
  }
  invisible(object)
}
