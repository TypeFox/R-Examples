summary.meta4diag = function(object,...){
  summarym4d = list()
  model.type = object$misc$model.type
  link = object$misc$link
  summarym4d$cpu.used = object$cpu.used
  summarym4d$summary.fixed = object[["summary.fixed"]]
  summarym4d$summary.hyperpar = object[["summary.hyperpar"]]
  summarym4d$summary.expected.accuracy = object$summary.expected.accuracy
  summarym4d$mlik = object$mlik[2,]
  summarym4d$var.type = c(rownames(object$summary.fixed),"var1", "var2", "rho")
  if(!object$misc$covariates.flag){
    summarym4d$correlation.expected.g = object$correlation.linear.comb
  }
  summarym4d$modality.name = object$misc$modality.name
  summarym4d$modality.level = object$misc$modality.level
  summarym4d$covariates.flag = object$misc$covariates.flag
  summarym4d$link = link
  summarym4d$model.type = model.type
  class(summarym4d) = "summary.meta4diag"
  return(summarym4d)
}

print.summary.meta4diag = function(x,...){
  cat('Time used: \n')
  print(x$cpu.used)
  cat('\n')
  cat('Fixed effects: \n')
  fixed = round(x[["summary.fixed"]],3)
  print(fixed)
  cat('\n')
  cat('Model hyperpar: \n')
  hyperpar = round(x[["summary.hyperpar"]],3)
  rownames(hyperpar)  = paste(rownames(x[["summary.hyperpar"]])," ",sep="")
  print(hyperpar)
  if(!x$covariates.flag){
    cat('\n')
    cat('-------------------')
    cat('\n')
    length = dim(fixed)[1]
    summarised.fixed = round(x$summary.expected.accuracy,3)
    print(summarised.fixed[1:length,])
  }
  if(!x$covariates.flag){
    cat('\n')
    cat('-------------------')
    cat('\n')
    if(is.null(x$modality.name)){
      cat(paste('Correlation between ',paste(rownames(x$summary.fixed),collapse=" and ")," is ",round(x$correlation.expected.g,4),".",sep=""))
    }else{
      paired.length = 0.5*dim(x$summary.fixed)[1]
      for(i in 1:paired.length){
        cat(paste('Correlation between ',paste(rownames(x$summary.fixed)[c(i, (i+paired.length))],collapse=" and ")," is ",round(x$correlation.expected.g[i],4),". \n",sep=""))
      }
    }
  }
  
  cat('\n')
  mlik = round(x$mlik,4)
  names(mlik) = ""
  cat(paste("Marginal log-likelihood: ",mlik,sep=""))
  cat('\n')
  cat('Variable names for marginal plotting: \n')
  cat("      ")
  cat(paste(x$var.type,collapse=", "))
  cat("\n")
}

