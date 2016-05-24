
residuals.HEGYtest <- function(object, ...) residuals(object$fit)

print.HEGYtest <- function(x, digits = max(3L, getOption("digits") - 2L), 
  dig.tst = max(1L, min(5L, digits - 1L)), print.fittedmodel = FALSE, ...)
{
  tab <- cbind(round(cbind(statistic=x$stat, "p-value"=x$pval), dig.tst), format(x$pvlabels))

  cat("\n")
  cat(strwrap(x$method, prefix = "\t"), sep = "\n")
  cat("\n")
  cat("data:  ", x$data.name, "\n\n", sep = "")

  if (print.fittedmodel)
  {
    x$fitted.model$series <- x$data.name
    cat("Fitted model\n")
    cat("------------\n")
    print(summary(x$fitted.model))
    cat("Test statistic\n")
    cat("--------------\n")
  }

  print(tab, quote = FALSE, right = TRUE, digits = digits, ...)
  cat("---\n")
  cat("Signif. codes:", attributes(x$pvlabels)$legend, "\n")

  cat("\n")
  cat(paste("Deterministic terms:", x$strdet, "\n"))
  cat(paste("Lag selection criterion and order: ", 
    x$lag.method, ", ", x$lag.order, "\n", sep = ""))
  cat(paste("P-values:", switch(x$type.pvalue, "RS" = "based on response surface regressions", 
    "bootstrap" = paste0("based on bootstrap (", x$bootstrap$nb, " replicates)"), 
    "raw" = "by interpolation in original tables"), "\n"))
}

summary.HEGYtest <- function(object, ...) { print.HEGYtest(object, print.fittedmodel = TRUE, ...) }
