print.bigssg <- 
  function(x,...){
    
    xnames <- names(x$xvars)
    nv <- length(xnames)
    if(nv==1L){
      cat("\nFamily:\n")
      cat(x$family,"\n")
      cat("\nSpline Type:\n")
      cat(paste(xnames,x$type),"\n")
      cat("\nFit Statistics:\n")
      print(x$info)
      cat("\nSmoothing Parameter:\n")
      cat(x$modelspec$lambda,"\n ")
    } else {
      cat("\nFamily:\n")
      cat(x$family,"\n")
      cat("\nSpline Types:\n")
      print(as.data.frame(x$type,row.names=""))
      cat("\nFit Statistics:\n")
      print(x$info)
      cat("\nAlgorithm Converged:\n")
      if(is.na(x$converge)){
        cat("Iterative update skipped (skip.iter=TRUE)")
      } else {
        cat(x$converged)
      }
      cat("\n ")
    } # end if(nv==1L)
    
  }