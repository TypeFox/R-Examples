print.bigspline <- 
  function(x,...){
    
    cat("\nSpline Type:\n")
    cat(x$type,"\n")
    cat("\nFit Statistics:\n")
    print(x$info)
    cat("\nSmoothing Parameter:\n")
    cat(x$lambda,"\n ")
    
  }