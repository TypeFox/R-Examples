summary.sparsenet=function(object,digits = max(3, getOption("digits") - 3),...){
  cat("\nCall: ", deparse(object$call), "\n\n")
  coeflist=object$coefficients
  rsq=object$rsq
  gnames=names(coeflist)
  dfmat=sapply(coeflist,"[[","df")
  list(Rsq=t(signif(rsq,digits)),Df=dfmat)
}
