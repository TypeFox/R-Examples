print.NSM3Ch5c <-function(x,...){
  
  if(x$method!="Exact"){
    cat("\n",x$method, " Approximation ")
      if(x$method=="Monte Carlo"){cat("(with ",x$n.mc, " Iterations) ")}
    cat("used: \n \n")
  }
  cat("Number of X values: ", x$m, "Number of Y values: ", x$n, "\n")
  
	if(x$method!="Asymptotic"){
	  if(!is.null(x$two.sided)){
      cat("For the given alpha=", x$alpha, ", the lower cutoff value is ",x$stat.name,"=",x$cutoff.L, ",\n", "with true alpha level=",round(x$true.alpha.L,4), "\n")
    }  
	 	cat("For the given alpha=", x$alpha, ", the upper cutoff value is ",x$stat.name, "=" ,x$cutoff.U, ",\n", "with true alpha level=",round(x$true.alpha.U,4), "\n")	
	}
	if(x$method=="Asymptotic"){
		if(!is.null(x$two.sided)){
			cat("For the given alpha=", x$alpha, ", the approximate lower cutoff value is ",x$stat.name,"=",x$cutoff.L, "\n")
		}		
		cat("For the given alpha=", x$alpha, ", the approximate upper cutoff value is ",x$stat.name, "=",x$cutoff.U, ",\n")
	}
	if(!is.null(x$extra)){
		cat(x$extra, "\n")
	}
}
