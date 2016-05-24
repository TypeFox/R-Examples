summary.LogConcDEAD <- function( object,...)
  { niter <- object$NumberOfEvaluations[1]
  
### What to do if got through getinfolcd;
  if( is.na( niter ) ) {
    warning( "Not computed by LogConcDEAD in this instance" )
    cat("\nLog MLE at observations: \n")
    print( object$logMLE);
  } else {
    ## other possibilities:
    if (niter==-1) errormessage <- "allocation error"
    else if (niter==-2) errormessage <- "improper space dimension"
    else if (niter==-3) errormessage <- "sigma returns an improper value"
    else if (niter==-4) errormessage <- "grad_sigma returns a zero or improper vector at the starting point"
    else if(niter==-7) errormessage <- "sigma is unbounded"
    else if (niter==-8) errormessage <- "gradient is zero, but stopping criteria are not fulfilled"
    else if (niter==-9) errormessage <- "iterations limit exceeded"
    else if (niter==-11) errormessage <- "Premature stop is possible"
    else if (niter==-12) errormessage <- "Result may not provide the true optimum"
    else if (niter==-13) errormessage <- "Result may be inaccurate in view of a point"
    else if (niter==-14) errormessage <- "Result may be inaccurate in view of a function value"
    
    if (niter > 0)
      {cat("\n Log MLE at observations: \n");
       print(object$logMLE);
       cat( "\n Number of Iterations: ",object$NumberOfEvaluations[1],"\n\n Number of Function Evaluations: ",object$NumberOfEvaluations[2],"\n")
     }
    else cat("SolvOpt error: ", errormessage, "\n")
  }
}


