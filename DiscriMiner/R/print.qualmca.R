#' @S3method print qualmca
print.qualmca <- 
function(x, ...)
{
  cat("\nMultiple Correspondence Analysis\n")
  cat(rep("-",43), sep="")
  cat("\n$values        ", "eigenvalues")
  cat("\n$coefficients  ", "coeffs of factorial axes")
  cat("\n$components    ", "factor coordinates\n")
  cat(rep("-",43), sep="")
  cat("\n\n$values\n")
  print(format(x$values, scientific=FALSE, digits=3),
        print.gap=2, quote=FALSE)
  cat("\n\n$coefficients\n")
  aux = min(5, ncol(x$coefficients))
  print(x$coefficients[,1:aux], print.gap=2)
  cat("...\n")
  invisible(x)
}
