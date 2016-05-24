print.NSM3Ch5p <-
function(x,...){
	if(!is.null(x$ties)&&x$ties){cat("Ties are present, so p-values are based on conditional null distribution. \n")}
  cat("Number of X values: ", x$m, "Number of Y values: ", x$n, "\n")
	cat(x$stat.name, " Statistic: ", round(x$obs.stat,4) , "\n")
	cat(x$method, " ")
	
  if(x$method=="Monte Carlo"){
		cat("(Using ", x$n.mc, "Iterations) ")
	}
  
	cat("upper-tail probability: ", ifelse(round(x$p.val,4)==0,x$p.val,round(x$p.val,4)), "\n")
	
  if(!is.null(x$two.sided)){
		cat(x$method, " ")
		if(x$method=="Monte Carlo"){
			cat("(Using ", x$n.mc, "Iterations) ")
		}
		cat("two-sided p-value: ", round(x$two.sided,4),"\n \n")
	}
	
  if(!is.null(x$extra)){
	  cat(x$extra, "\n")
	}
  
}
