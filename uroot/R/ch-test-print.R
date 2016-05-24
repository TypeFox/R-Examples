
print.CHtest <- function(x, digits = max(3L, getOption("digits") - 2L), 
  dig.tst = max(1L, min(5L, digits - 1L)), print.fittedmodel = FALSE, ...)
{
  tab <- cbind(round(cbind(statistic=x$stat, pvalue=x$pval), dig.tst), format(x$pvlabels))

  cat("\n")
  cat(strwrap(x$method, prefix = "\t"), sep = "\n")
  cat("\n")
  cat("data:  ", x$data.name, "\n\n", sep = "")

  if (print.fittedmodel)
  {
    x$fitted.model$series <- x$data.name
    cat("Fitted model\n")
    cat("------------\n")
    if (x$method.fit == "qr") print(summary(x$fitted.model)) else print(x$fitted.model)
    if (x$method.fit == "qr") cat("Test statistic\n") else cat("\nTest statistic\n")
    cat("--------------\n")
  }

  print(tab, quote = FALSE, right = TRUE, digits = digits, ...)
  cat("---\n")
  cat("Signif. codes:", attributes(x$pvlabels)$legend, "\n")

  if (!x$isNullxreg) {
    xregnms <- names(coef(x$fitted.model))
    xregnms <- xregnms[-seq_len(tail(grep("ypi.Ypi", xregnms),1))]
    xregnms <- paste(xregnms,collapse=", ")
  } else xregnms <- NULL

  cat("\n")
  cat(paste("Test type:",
    switch(x$type, "dummy" = "seasonal dummies", "seasonal cycles"), "\n"))
  cat(paste("NW covariance matrix lag order:", x$NW.order, "\n"))
  cat(paste("First order lag:", c("no", "yes")[as.numeric(x$lag1)+1], "\n"))
  cat(paste("Other regressors:", c("yes, ", "no")[as.numeric(x$isNullxreg)+1], xregnms, "\n"))
  cat(paste("P-values:", switch(x$type.pvalue, "RS" = "based on response surface regressions", 
    "raw" = "interpolation in original tables"), "\n"))
}

summary.CHtest <- function(object, ...) { print.CHtest(object, print.fittedmodel = TRUE, ...) }
