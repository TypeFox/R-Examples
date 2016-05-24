#' @S3method print disqual
print.disqual <- 
function(x, ...)
{
  cat("\nDiscriminant Analysis on Qualitative Variables\n")
  cat(rep("-",46), sep="")
  cat("\n$raw_coefs       ", "raw coeffcients")
  cat("\n$norm_coefs      ", "normalized coefficients")
  cat("\n$confusion       ", "confusion matrix")
  cat("\n$scores          ", "scores")
  cat("\n$classification  ", "assigned class")
  cat("\n$error_rate      ", "error rate\n")
  cat(rep("-",46), sep="")
  cat("\n\n$raw_coefs\n")
  print(format(x$raw_coefs, scientific=FALSE, digits=3),
        print.gap=2, quote=FALSE)
  cat("\n\n$norm_coefs\n")
  print(format(x$norm_coefs, scientific=FALSE, digits=3),
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


