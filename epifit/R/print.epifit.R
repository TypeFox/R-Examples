##' Print function for epifit object.
##'
##' This function print result of function \code{\link{epifit}}
##'
##' @param x Object of class \code{epifit}.
##' @param digits a non-null value for digits specifies the minimum number of significant digits to be printed in values. The default,  uses \code{max(\link[base:options]{getOption}}(digits - 4, 3)).
##' @param ... Further arguments passed to or from other methods.
##' @export
print.epifit <- function(x, digits=max(options()$digits -4, 3), ...){
  if(x$convergence==3){
    warning("Did not converge in nlm: last global step failed to locate a point (local minimum or 'stol is too slamm).\n")
  } else if(x$convergence==4){
    warning("Iteration limit exeeded in nlm.\n")
  } else if(x$convergence==5){
    warning("maximum step size exceeded five consecutive times in nlm (stepmax is too small or the function is unbounded below, becomes asymptotic to a finite value from above in some direction).\n")
  } else if(x$convergence==6){
    warning("Degeneracy of the Nelder-Mead simplex in optim.\n")
  } else if(x$convergence==7){ # not implemented
    warning("A warning from the L-BFGS-B method in optim.\n")
  } else if(x$convergence==8){ # not implemented
    warning("an error from the L-BFGS-B method in optim.\n")
  }
  cat("Call:\n")
  dput(x$call)
  cat("\n")
  se <- sqrt(diag(x$var))
  z <- x$coefficient/se
  p <- 2*(1 - pnorm(abs(z)))
  #logtest <- -2*(x$loglik[1] - x$loglik[2])
  df <- sum(!is.na(length(x$coefficient)))
  tbl <- cbind(x$coefficient, se, z, p)
  rownames(tbl) <- x$parameters
  colnames(tbl) <- c("coef", "se(coef)", "z", "p")
  print(tbl, digits = digits)
  #cat("\n")
  #cat("Likelihood ratio test=", format(round(logtest, 2)), "  on ",
  #    df, " df,", " p=", format(1 - pchisq(logtest, df)),"\n",  sep="")
  invisible(x)
}
