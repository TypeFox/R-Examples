print.NSM3Ch7c <-
  function(x,...){
    if(x$method!="Exact"){
      cat("\n",x$method, " Approximation ")
      if(x$method=="Monte Carlo"){cat("(with ",x$n.mc, " Iterations) ")}
      cat("used: \n \n")
    }
    
    if(is.null(x$trt)){
      cat("Number of blocks: n=", x$n,"\n")
      cat("Number of treatments: k=", x$k, "\n")
    }
    if(!is.null(x$ss)){
      cat("Number of treatments per block: s=", x$ss,"\n")
    }
    if(!is.null(x$pp)){
      cat("Number of observations per treatment: p=", x$pp,"\n")
    }
    if(!is.null(x$lambda)){
      cat("Number of times each pair of treatments occurs together within a block: lambda=", x$lambda,"\n")
    }
    
    if(!is.null(x$trt)){
      cat("Control group size: ", x$trt, "Treatment group size(s): ", x$n, "\n")
    }
    
    if(x$method!="Asymptotic"){  
      cat("For the given alpha=", x$alpha, ", the upper cutoff value is ",x$stat.name, "=" ,x$cutoff.U, ", with true alpha level=",round(x$true.alpha.U,4), "\n")  
    }
    if(x$method=="Asymptotic"){
      cat("For the given alpha=", x$alpha, ", the approximate upper cutoff value is ",x$stat.name, "=",x$cutoff.U, ",\n")
    }
    if(!is.null(x$extra)){
      cat(x$extra, "\n")
    }
  }