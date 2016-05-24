print.summary.ecoML <- function(x, digits=max(3,
                                     getOption("digits")-3), ...) {

  cat("\nCall: ", paste(deparse(x$call), sep="\n", collapse="\n"))
  cat("\n")
  if (!is.null(x$param.table)) {
    cat("\n*** Parameter Estimates ***\n")
    if (x$fix.rho)
      cat("\nOriginal Model Parameters (rho is fixed at ", x$rho, "):\n", sep="")   
    else
      cat("\nOriginal Model Parameters:\n")
    print(x$param.table, digits=digits, na.print="NA",...)
  }

  cat("\n*** Insample Predictions ***\n")
  cat("\nUnweighted:\n")
  print(x$agg.table, digits=digits, na.print="NA",...)
  
  if (!is.null(x$agg.wtable)) {
  cat("\nWeighted:\n")
  print(x$agg.wtable, digits=digits, na.print="NA",...)
  }
  if (!is.null(x$W.table)) {
    cat("\n\nUnit-level Estimates of W:\n")
    print(x$W.table, digits=digits, na.print="NA",...)
  }

  cat("\n\nLog-likelihood:", x$loglik)
  cat("\nNumber of Observations:", x$n.obs)
  cat("\nNumber of EM iterations:", x$iters.em)
  if (x$sem)
    cat("\nNumber of SEM iterations:", x$iters.sem)
  cat("\nConvergence threshold for EM:", x$epsilon)
  
  cat("\n\n")
  invisible(x)
}
