print.icweib <-
function(x, digits=3, ...) {
  cat("Total observations used:", x$ns[1])
  cat(". Model Convergence:", x$q$convergence==0, "\n")
  if (x$ns[4]>0) cat(x$ns[4], " observations are deleted due to inappropriate times or missing data\n",
                  "NOTE: 0/Inf are used for left/right censoring instead of NA\n", sep="")
  if (x$ns[3] > 0) {
    cat("\nCoefficients: \n")
    print(x$coef, digits=digits, ...)  	
  }
  cat("\nWeibull parameters - gamma(shape), lambda(scale): \n")
  print(x$weib, digits=digits, row.names=F, ...)
  if (x$ns[2]>1) {	
  cat("\nTest of proportional hazards for strata (H0: all strata's shape parameters are equal):\n")
  print(x$stratatest, digits=digits, row.names=F, ...)
  }
  cat("\n")
  cat("Loglik(model)= ", x$loglik[1])
  if (x$ns[2] > 1) cat("   Loglik(reduced)= ", x$loglik[2])
  cat("\nLoglik(null)= ", x$loglik[3])
  if (x$ns[3] > 0) {
  	df <- x$ns[3]
  	likratio <- 2*(x$loglik[1] - x$loglik[3])
  	plik <- 1 - pchisq(likratio, df)
  	cat("  Chisq=", likratio, "  df=", df, " p.value=", plik)
  }
}
