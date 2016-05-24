print.bgeva <- function(x,...){

  cat("\nFamily: BGEVA \nEquation: ")
  print(x$gam.fit$formula)

  cat("\n")

  cat("n = ",x$n,"  tau = ",format(x$tau, digits=3),"  total edf = ",format(x$t.edf, digits=3),"\n\n")

invisible(x)

}
