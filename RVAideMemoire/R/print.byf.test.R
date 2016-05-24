print.byf.test <- function(x,...) {
  cat("\n")
  cat(strwrap(x$method,prefix="\t"),sep="\n")
  cat("\n")
  cat("data: ",x$data.name,"\n\n")
  printCoefmat(x$tab,digits=5,P.values=TRUE,has.Pvalue=TRUE)
  cat("\n")
}
