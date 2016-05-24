#' @S3method print linda
print.linda <- 
function(x, ...)
{
  cat("\nLinear Discriminant Analysis\n")
  cat(rep("-",43), sep="")
  cat("\n$functions       ", "discrimination functions")
  cat("\n$confusion       ", "confusion matrix")
  cat("\n$scores          ", "discriminant scores")
  cat("\n$classification  ", "assigned class")
  cat("\n$error_rate      ", "error rate\n")
  cat(rep("-",43), sep="")
  cat("\n\n$functions\n")
  print(format(x$functions, scientific=FALSE, digits=3),
        print.gap=2, quote=FALSE)
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
