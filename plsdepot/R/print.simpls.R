#' @S3method print simpls
print.simpls <-
function(x,...)
{       
  cat("\nSIMPLS\n")
  cat(rep("-",41), sep="")
  cat("\n$x.scores   ", "X-scores (T-components)")
  cat("\n$x.wgs      ", "X-weights")
  cat("\n$y.wgs      ", "Y-weights")
  cat("\n$cor.xt     ", "X,T correlations")
  cat("\n$cor.yt     ", "Y,T correlations")
  cat("\n$R2X        ", "explained variance of X by T")
  cat("\n$R2Y        ", "explained variance of Y by T\n")
  cat(rep("-",41), sep="")
  cat("\n\n")
  invisible(x)
}