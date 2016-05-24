#' @S3method print quada
print.quada <- 
function(x, ...)
{
  cat("\nQuadratic Discriminant Analysis\n")
  cat(rep("-",37), sep="")
  cat("\n$confusion       ", "confusion matrix")
  cat("\n$scores          ", "discriminant values")
  cat("\n$classification  ", "assigned class")
  cat("\n$error_rate      ", "error rate\n")
  cat(rep("-",37), sep="")
  cat("\n\n$confusion\n")
  print(x$confusion, print.gap=2)
  cat("\n\n$error_rate\n")
  print(x$error_rate)
  cat("\n\n$scores\n")
  print(head(x$scores), print.gap=2)
  cat("...\n")
  cat("\n$classification\n")
  print(head(x$classification))
  cat("...\n")
  invisible(x)
}
