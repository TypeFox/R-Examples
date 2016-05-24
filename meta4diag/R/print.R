print.meta4diag = function(x,...){
  cat('Time used: \n')
  print(x$cpu.used)
  cat("\n")
  cat(paste('Model:Binomial-Normal Bivariate Model for ', paste(x$names.fitted,collapse=" & "),". \n",sep=""))
  cat(paste('Data contains ', dim(x$data)[1], " primary studies. \n",sep=""))
  cat("\n")
  if(x$misc$modality.flag){
    cat(paste("Data has Modality variable with level ", x$misc$modality.level, ". \n",sep=""))
    if(x$misc$covariates.flag){
      cat(paste("Data has Covariates variable with name ", x$misc$covariates.name, ". \n",sep=""))
    } else{
      cat("Covariates not contained. \n")
    }  
  }else{ # no modality
    cat("Data has no Modality variable. \n")
    if(x$misc$covariates.flag){
      cat(paste("Data has Covariates variable with name ", x$misc$covariates.name, ". \n",sep=""))
    } else{
      cat("Covariates not contained. \n")
    }
  }
  cat("\n")
  cat(paste("Model using link function ",x$misc$link,".\n",sep=""))
  cat("\n")
  cat("Marginals can be plotted with setting variable names to ") 
  cat("\n")
  cat(paste(paste(rownames(x$summary.fixed),collapse=", "), ", var1, var2 and rho. \n",sep=""))
  
  return(invisible())
}