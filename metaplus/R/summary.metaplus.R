summary.metaplus <- function(object, ...) {
  if (!inherits(object, "metaplus"))
    stop("Use only with 'metaplus' objects.\n")

  if (object$justfit) stop("Cannot use with objects fitted with justfit=TRUE")
  
  fitstats <- list(logLik=logLik(object),AIC=AIC(object),BIC=BIC(object))
  out <- list(results=object$results,fitstats=fitstats)
  class(out) <- "summary.metaplus"
  out
}

print.summary.metaplus <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!inherits(x, "summary.metaplus"))
    stop("Use only with 'metaplus' objects.\n")
  #  print.default(object$results,na.print="")
   printCoefmat(x$results,digits=digits, signif.stars=FALSE,cs.ind=1:3,tst.ind=integer(),P.values=TRUE,has.Pvalue=TRUE,
               eps.Pvalue=1.0e-4,na.print="",...)
  cat("\n")
  print(data.frame(logLik = x$fitstats$logLik,AIC = x$fitstats$AIC, BIC = x$fitstats$BIC,
                     row.names = " "))
  invisible(x)
}