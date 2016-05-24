print.Mort1Dsmooth <-
function(x,
                               digits=max(3,
                                 getOption("digits")-3),
                               ...){
  if(!is.null(cl <- x$call)){
    cat("Call:\n")
    dput(cl, control=NULL)
  }
  cat("\nNumber of Observations              :", x$n, "\n")
  cat("Effective dimension                 :", format(signif(x$df, 4)), "\n")
  cat("(Selected) smoothing parameter      :", format(signif(x$lambda, 5)), "\n")
  cat("Bayesian Information Criterion (BIC):", signif(x$bic, 5))
  cat("\n")
  invisible(x)
}
