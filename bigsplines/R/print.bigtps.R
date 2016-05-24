print.bigtps <- 
  function(x,...){
    
    cat("\nSpline Type:\n")
    cat("tps","\n")
    cat("\nFit Statistics:\n")
    print(x$info)
    cat("\nSmoothing Parameter:\n")
    cat(x$lambda,"\n ")
    
  }