# last modified 12 Dec 04 by J. Fox

"print.hetcor" <-
function(x, digits = max(3, getOption("digits") - 3), ...){
  R <- signif(x$correlations, digits=digits)
  R[upper.tri(R)] <- x$type[upper.tri(R)]
  R <- as.data.frame(R)
  if (x$ML) cat("\nMaximum-Likelihood Estimates\n")
  else cat("\nTwo-Step Estimates\n")
  cat("\nCorrelations/Type of Correlation:\n")
  print(R)
  if (!is.null(x$std.errors)){
    SE <- signif(x$std.errors, digits)
    diag(SE) <- ""
    if (x$NA.method == "complete.obs"){
      SE[upper.tri(SE)] <- ""
      cat("\nStandard Errors:\n")
      SE <- as.data.frame(SE)
      print(SE[,-ncol(SE)])
      cat(paste("\nn =", x$n, "\n"))
      }
    else {
      SE[upper.tri(SE)] <- x$n[upper.tri(SE)]
      diag(SE) <- diag(x$n)
      SE <- as.data.frame(SE)
      cat("\nStandard Errors/Numbers of Observations:\n")
      print(SE)
      }
    if (!all(is.na(x$tests[lower.tri(x$tests)]))){
      Test <- signif(x$tests, digits)
      Test[upper.tri(Test)] <- ""
      diag(Test) <- ""
      Test <- as.data.frame(Test)
      cat("\nP-values for Tests of Bivariate Normality:\n")
      print(Test[,-ncol(Test)])
      }
    }
  invisible(x)
  }
