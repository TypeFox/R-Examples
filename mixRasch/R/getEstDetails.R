`getEstDetails` <- function(raschResult, camelCase=TRUE){  
    
	if(! is.array(raschResult$class)){
  	  n.c <- 1
	} else{
	  n.c <- ncol(raschResult$class)
    }	
	
	out <- list(model = raschResult$model, n.c = n.c, iter = raschResult$iter, max.change=raschResult$max.change, converge.flag=raschResult$converge.flag, run.time=raschResult$run.time)	
    
    if(camelCase){
      names(out) <- c("model", "nC", "iter", "maxChange", "convergeFlag", "runTime")
    }		  	  
    out
  }