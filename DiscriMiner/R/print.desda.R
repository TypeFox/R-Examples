#' @S3method print desda
print.desda <-
function(x, ...)
{
  cat("\nDescriptive Discriminant Analysis\n")
  cat(rep("-",33), sep="")
  cat("\n$power     ", "discriminant power")
  cat("\n$values    ", "table of eigenvalues")
  cat("\n$discrivar ", "discriminant variables")
  cat("\n$discor    ", "correlations")
  cat("\n$scores    ", "discriminant scores\n")
  cat(rep("-",33), sep="")
  cat("\n\n$power\n")
  print(format(x$power, scientific=FALSE, digits=4),
        print.gap=2, quote=FALSE)
  cat("\n\n$values\n")
  print(format(x$values, scientific=FALSE, digits=3),
        print.gap=2, quote=FALSE)
  cat("\n\n$discrivar\n")
  print(x$discrivar, print.gap=2, digits=4)
  cat("\n\n$discor\n")
  print(head(x$discor), print.gap=2, digits=4)
  cat("\n\n$scores\n")
  print(head(x$scores), print.gap=2, digits=4)
  cat("...\n\n")
  invisible(x)
}
