print.summary.npregress <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  resid <- x$residuals
  df <- x$df
  nam <- c("Min", "1Q", "Median", "3Q", "Max")
  rq <-structure(quantile(resid), names = nam)
  cat("Residuals:\n")
  print(rq, digits = digits, ...)
  cat("Residual standard error:", format(signif(x$Std.Error, 
        digits)), "on", format(round(x$Resid.Df,1)), "degrees of freedom\n")
  if (x$criteria!="user") {
    crite <- structure(x$criteria, names = names(x$criteria))
    print(crite, digits = digits, ...)
  }
  kernelsmooth <- c("gaussian","Epanechnikov","uniform","quartic")
  resu <- match(x$kernel,c("g","e","u","q"))
  cat("Kernel:",kernelsmooth[resu],"(with",format(signif(x$Df, 
        digits)),"df)\n")
  cat("\nBandwidth:", format(signif(x$bandwidth,digits)), "chosen by", x$crit4bw, "\n")
  invisible(x)
}
