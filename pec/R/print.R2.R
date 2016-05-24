
##' @export
print.R2 <- function(x,...){
  cat("\nTime-dependent explained variation:\n\n 1- Brier(model)/Brier(reference)\n\nReference model: ",attr(x,"reference"),"\n\n")
  print.listof(x,...)
}
