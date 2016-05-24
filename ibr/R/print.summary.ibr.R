print.summary.ibr <- function(x, displaybw=FALSE, digits = max(3, getOption("digits") - 3), ...) {
  resid <- x$residuals
  df <- x$df
  nam <- c("Min", "1Q", "Median", "3Q", "Max")
  rq <-structure(quantile(resid), names = nam)
  cat("Residuals:\n")
  print(rq, digits = digits, ...)
  cat("Residual standard error:", format(signif(x$Std.Error, 
        digits)), "on", format(round(x$Resid.Df,1)), "degrees of freedom\n")
  cat("\nInitial df:", format(round(x$Initial.Df,2)), "; Final df:", format(round(x$Final.Df,2)), "\n")
  if (!is.null(x$criteria)) {
      crite <-structure(x$criteria, names = names(x$criteria))
      if (length(crite)>1) {
          cat("\nCriteria:\n")
      }
      print(crite, digits = digits, ...)
  if ((length(crite)>1)&&(!all(is.na(x$iterations)))) {
    cat("Iterations:\n")
  crite2 <- structure(x$iterations,names = names(x$iterations))
    print(crite2, digits = digits, ...)
}
  }
  cat("\nNumber of iterations:", x$iter, "chosen by", x$crit4iter, "\n")
  kernelsmooth <- c("gaussian kernel","Epanechnikov kernel","uniform kernel","quartic kernel")
  resu <- match(x$kernel,c("g","e","u","q"))
  if (x$smoother=="k") {
    cat("Base smoother:",kernelsmooth[resu],"(with",format(round(x$Initial.Df,2)),"df)\n")
    if (displaybw) {
      cat("\nBandwith\n")
      bw <-structure(x$bandwidth, names = names(x$bandwidth))
      print(bw, digits = digits, ...)
    }
  }
  if (x$smoother=="tps") {
    cat("Base smoother: Thin plate spline of order",x$m,"(with",format(round(x$Initial.Df,2)),"df)\n")
    if (displaybw) {
      cat("\nBandwith\n")
      bw <-structure(x$bandwidth, names = "lambda")
      print(bw, digits = digits, ...)
    }
  }
  if (x$smoother=="ds") {
    cat("Base smoother: Duchon spline with derivative order m=",x$m,"\n  and weight exponent s=",x$s,"(with",format(round(x$Initial.Df,2)),"df)\n")
    if (displaybw) {
      cat("\nBandwith\n")
      bw <-structure(x$bandwidth, names = "lambda")
      print(bw, digits = digits, ...)
    }
  }
  invisible(x)
}
