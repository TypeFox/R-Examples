print.isoph=function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nEstimated monotone function, hazard ratio and level set:\n")
  print(x$est)
  if(!is.na(x$exp.beta)){
    a1.0=paste("\n exp(beta): ",x$exp.beta,"\n",sep='')
    cat(a1.0)
  }

  a1.1=paste("\nNumber of events/subjects was ",x$nevent,"/",x$n,".",sep='')
  a1.2=paste("\nNumber of distinct covariates associated with observed events was ",x$njump,".",sep='')
  a1.3=paste("\nShape restriction was monotone ",x$shape,".",sep='')
  cat(a1.1)
  cat(a1.2)
  cat(a1.3)
  
  if(x$conv=="converged"){
    a2=paste("\n\nAlgorithm was ",x$conv,".\n",sep='')
  }else if(x$conv=="not converged"){
    a2=paste("\n\nAlgorithm was ",x$conv,". Results were questinable.\n",sep='')
  }
  cat(a2)
}


