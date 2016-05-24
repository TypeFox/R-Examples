print.uniah=function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nEstimated shape restricted hazard function and level set:\n")
  print(x$est)
  if(!is.na(x$beta)){
    a1.0=paste("\n beta: ",x$beta,"\n",sep='')
    cat(a1.0)
  }
  
  if(is.null(x$lf)){ #known
    a1.1=paste("\nMode was set to ",x$M,sep='')
  }else{             #unknown mode
    a1.1=paste("\nEstimated Mode: ",x$M,sep='')
  }
  
  a1.2=paste("\nNumber of events/subjects was ",x$nevent,"/",x$n,".",sep='')
  a1.3=paste("\nNumber of distinct covariates associated with observed events was ",x$njump,".",sep='')
  a1.4=paste("\nShape restriction was ",x$shape,".",sep='')

  cat(a1.1)
  cat(a1.2)
  cat(a1.3)
  cat(a1.4)  
  
  if(x$conv=="converged"){
    a2=paste("\n\nAlgorithm was ",x$conv, ".\n",sep='')
  }else{
    a2=paste("\n\nAlgorithm was ",x$conv,". Results were questinable.\n",sep='')
  }
  cat(a2)
  
}



