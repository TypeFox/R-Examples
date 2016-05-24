#' @S3method print nipals
print.nipals <-
function(x, ...)
{
  cat("\nNIPALS algorithm\n")
  cat(rep("-",34), sep="")
  cat("\n$values    ", "eigenvalues")
  cat("\n$scores    ", "scores (T-components)")
  cat("\n$loadings  ", "loadings")
  cat("\n$cor.xt    ", "X,T correlations")
  cat("\n$disto     ", "distance to origin")
  cat("\n$contrib   ", "contribution of rows")
  cat("\n$cos       ", "squared cosinus")
  cat("\n$dmod      ", "distance to the model\n")
  cat(rep("-",34), sep="")
  cat("\n\n")
  invisible(x)
}

