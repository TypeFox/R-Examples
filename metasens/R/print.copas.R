print.copas <- function(x, sign.rsb=x$sign.rsb,
                        backtransf=x$backtransf,
                        digits=max(3, .Options$digits - 3),
                        ...){
  
  meta:::chkclass(x, "copas")
  
  
  cl <- class(x)[1]
  addargs <- names(list(...))
  ##
  fun <- "print.copas"
  ##
  meta:::warnarg("logscale", addargs, fun, otherarg="backtransf")
  ##
  if (is.null(backtransf))
    if (!is.null(list(...)[["logscale"]]))
      backtransf <- !list(...)[["logscale"]]
    else
      backtransf <- TRUE
  ##
  if (is.null(sign.rsb))
    sign.rsb <- 0.1
  else
    meta:::chklevel(sign.rsb)
  
  
  meta:::crtitle(x)

  cat("Copas selection model analysis\n\n")
  
  
  res <- cbind(c("range of gamma0: ", "range of gamma1: "),
               format(c(round(x$gamma0.range[1], 2),
                        round(x$gamma1.range[1], 2))),
               ##               c(" - ", " - "),
               format(c(round(x$gamma0.range[2], 2),
                        round(x$gamma1.range[2], 2))))
  ##
  dimnames(res) <- list(rep("", dim(res)[1]),
                        c("", "min", "max"))
  ##
  prmatrix(res, quote=FALSE, right=TRUE)
  
  
  cat("\nLargest standard error (SE):", max(round(x$seTE, 4)), "\n\n")
  ##
  cat("Range of probability publishing trial with largest SE:\n")
  ##
  res <- matrix(format(round(range(pnorm(x$gamma0+x$gamma1/max(x$seTE))),
                             3)), nrow=1)
  ##
  dimnames(res) <- list(rep("", dim(res)[1]), c("min", "max"))
  ##
  prmatrix(res, quote=FALSE, right=TRUE)

  cat("\n\nCalculation of orthogonal line:\n\n")
  ##
  res <- as.matrix(data.frame(x$regr)[,c("levels", "nobs",
                                         "adj.r.squareds",
                                         "slopes", "se.slopes")])
  dimnames(res) <- list(rep("", dim(res)[1]),
                        c("level", "nobs",
                          "adj.r.square",
                          "slope", "se.slope"))
  prmatrix(res, quote=FALSE, right=TRUE)
  
  cat("\n\n")
  print(summary(x, sign.rsb=sign.rsb),
        digits=digits, header=FALSE, backtransf=backtransf)
  
  invisible(NULL)
}
