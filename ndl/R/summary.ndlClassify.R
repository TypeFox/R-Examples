summary.ndlClassify <- function(object, ...)
{
  statistics <- ndlStatistics(object, ...)

  weights <- object$weightMatrix

  sumry <- list(call=object$call, formula=object$formula, weights=weights, statistics=statistics, ...)

  class(sumry) <- "summary.ndlClassify"
  return(sumry)
}

print.summary.ndlClassify <- function(x, digits=max(3,getOption("digits")-3), max.print=10, ...)
{
  if(is.na(max.print))
    max.print= NROW(x$weights)
#  if(!is.null(x$digits) & is.numeric(x$digits))
#    digits=x$digits
#  if(!is.null(x$max.print) & is.numeric(x$max.print))
#    max.print=x$max.print

  cat("\nCall:\n")
  print(x$call)
  cat("\nFormula:\n")
  print(x$formula)
  cat("\nWeights:\n")
  tabl <- x$weights[1:min(nrow(x$weights),max.print),]
  class(tabl) <- "table"
  print(tabl, digits=digits)
  if(nrow(x$weights)>max.print)
    cat(paste("... [ omitted ",nrow(x$weights)-max.print," rows ] ...\n",sep=""))
  cat("\n")
  statistics <- x$statistics
  deviances <- format(c(signif(statistics$deviance.null,digits),signif(statistics$deviance.model,digits)))
  DFs <- format(c(statistics$df.null,statistics$df.model))
  cat(c("Null deviance:             ", deviances[1], " on ", DFs[1], " degrees of freedom\n"))
  cat(c("Residual (model) deviance: ", deviances[2], " on ", DFs[2], " degrees of freedom\n"))
  cat(format(c("\nR2.likelihood: ", signif(statistics$R2.likelihood,digits), "\nAIC: ", signif(statistics$AIC.model,digits), "\nBIC: ", signif(statistics$BIC.model,digits))))
  cat("\n\n")

  invisible(x)
}
