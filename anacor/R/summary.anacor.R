`summary.anacor` <-
function(object, ...)
{
# object ... object of class "anacor"
  cat("\nz-test for singular values:\n")
  z.values <- object$singular.values/object$se.singular.values
  pvalues <- (1-pnorm(abs(z.values)))
  stestmat <- data.frame(round(object$singular.values,4),round(object$se.singular.values,4),round(pvalues,4))
  colnames(stestmat) <- c("Singular Values","Asymptotical SE","p-value")
  rownames(stestmat) <- paste("D", 1:object$ndim, sep = "")
  print(stestmat)
  cat("\nRow scores:\n")
  print(round(object$row.scores,5))
  cat("\nColumn scores:\n")
  print(round(object$col.scores,5))
  cat("\n")
}

