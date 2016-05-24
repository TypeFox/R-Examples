print.NSM3Ch6MCc <-
  function(x,...){
    if(x$method!="Exact"){
      cat("\n",x$trtethod, " Approximation ")
      if(x$method=="Monte Carlo"){cat("(with ",x$n.mc, " Iterations) ")}
      cat("used: \n \n")
    }
    
    if(is.null(x$trt)){
      cat("Group sizes: ", x$n, "\n")
    }
    
    if(!is.null(x$trt)){
      cat("Control group size: ", x$trt, "Treatment group size(s): ", x$n, "\n")
    }
    
    if(x$method!="Asymptotic"){  
      cat("For the given experimentwise alpha=", x$alpha, ", the upper cutoff value is",x$stat.name, "=" ,x$cutoff.U, ",", "with true experimentwise alpha level=",round(x$true.alpha.U,4), "\n")  
    }
    if(x$method=="Asymptotic"){
      cat("For the given experimentwise alpha=", x$alpha, ", the approximate upper cutoff value is",x$stat.name, "=",x$cutoff.U, ",\n")
    }
    if(!is.null(x$extra)){
      cat(x$extra, "\n")
    }
}
