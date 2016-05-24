

print.hzar.fitRequest <- function(x,noLLFB=TRUE,...){
  cat("Fit request object:\n")
  ## Model Param
  cat("\nModel Parameters:\n")
  mP=x$modelParam
  print(data.frame(row.names=names(mP$init),
                   init=as.numeric(mP$init),
                   tune=as.numeric(mP$tune),
                   lower=as.numeric(mP$lower),
                   upper=as.numeric(mP$upper)),...)
  cat("\nFixed Parameters:\n")
  cat(paste(names(mP$fixed),"=",as.numeric(mP$fixed),collapse="\n"))
  cat("\n")
  cat("\nMCMC properties:\n")
  cP=x$mcmcParam
  for(iter in names(cP)){
    if(iter == "seed" && length(cP[[iter]])>1){
      cat(paste("  seed base =","[",paste(cP[[iter]][[1]],collapse=", "),"]\n"))
      cat(paste("    channel =",paste(cP[[iter]][[2]],collapse=", "),"\n"))
    } else {
      cat(paste(sprintf("%11s",iter),"=",paste(cP[[iter]],collapse=", "),"\n"))
    }
  }
  
  if(!noLLFB){
    cat("\nCompiled LL expression:\n")
    llB <- body(x$llFunc)
    cat(paste(deparse(llB),collapse="\n"))
  } else {
    cat("\nLL expression hidden\n")
  }
  if(!is.null(x$cM)){
    if(nrow(x$cM)<5){
      cat("\nCovariance Matrix:\n")
      print(x$cM)
    } else {
      cat("\nCovariance Matrix Diagonal:\n")
      print(diag(x$cM))
    }
  }
  if(!is.null(x$mcmcRaw)&&nrow(x$mcmcRaw)>1){
##     if(ncol(x$mcmcRaw)<5){
    cat("\nMCMC trace summary:\n")
    print(summary(x$mcmcRaw))
  }
  if(attr(x,"fit.run")){
    cat("\nFit request has been run ")
    if(attr(x,"fit.success"))
      cat("successfully.\n")
    else
      cat("unsuccessfully.\n")
  }
      
}
