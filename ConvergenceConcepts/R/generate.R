# We generate the sample paths
##############################

# randomgen: it is the generator of X_n-X values, or only of X_n values in case we want to study convergence in law

# argsgen should be given as a list

generate <- function(randomgen,nmax=1000,M=500,argsgen=NULL) {
  
  if(is.null(argsgen)) {
    
    data <- matrix(replicate(M,randomgen(nmax)),nrow=M,ncol=nmax,byrow=TRUE)
    
  }
  
  else {
    res <- rep("",length(argsgen))

    for (i in 1:length(argsgen)) res[i] <- paste(names(argsgen)[i],"=",argsgen[i],sep="")
    
    argsgen <- paste(res,collapse=",")
    
    eval(parse(text=paste("data <- matrix(replicate(M,randomgen(nmax,",argsgen,")),nrow=M,ncol=nmax,byrow=TRUE)",sep="")))
    
  }
  
  return(list(data=data))

}


