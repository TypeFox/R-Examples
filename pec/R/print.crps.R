##' @export
print.crps <- function(x,digits=3,...){
  cat("\nIntegrated Brier score (crps):\n\n")
  if (is.list(x))
    print.listof(x,digits=digits,...)
  else
    prmatrix(round(x,digits=3),...)
}
