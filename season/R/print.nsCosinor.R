## print.nsCosinor.R
## Prints basic results from nsCosinor

print.nsCosinor<-function(x, ...){

  ## Checks
  if (class(x)!="nsCosinor"){
    stop("Object must be of class 'nsCosinor'")} 

  ## Statistics ###
  cat("Non-stationary cosinor\n\n")
  cat("Call:\n")
  print(x$call)
  cat("\nNumber of MCMC samples = ",x$call$niters-x$call$burnin+1,
      "\n\n",sep="")
  cat("Length of time series = ",x$n,"\n",sep="")
  cat("\nResidual statistics\n",sep="")
  print(summary(x$residuals), ...)
} # end of function
