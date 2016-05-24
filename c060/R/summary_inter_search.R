
###########################################################################################################
summary.intsearch<-function(object,digits = max(3, getOption("digits") - 3), verbose=TRUE, first.n=5, ...){
  fit <- object
  alphas <- fit$Xtrain[,1]
  lambdas <- unlist(sapply(sapply(fit$model, "[", "model"), "[", "lambda"))
  deviances <- fit$Ytrain
  # round problems!!!! take first from the fit 
  # number of selected features in the models; dfs
  tmp.models<-sapply(sapply(sapply(fit$model, "[", "model"), "[", "cvreg"), "[", "glmnet.fit")
  
  n.features<-mapply( function(List, lam) List$df[which(List$lambda %in% lam)], tmp.models, lambdas)
  
  # optimal models
  #print("chose the model with min num of FS ")      
  opt.models <- sapply(fit$model.list, "[", "model") [fit$Ytrain == fit$fmin ]
  
  opt.alpha <- opt.models[[1]]$alpha
  opt.lambda <- opt.models[[1]]$lambda
  opt.error <- fit$fmin 
  
  out <- list(info=data.frame(alpha=alphas,lambda=lambdas,deviance=deviances,n.features=n.features),
              opt.alpha=opt.alpha, opt.lambda=opt.lambda, opt.error=opt.error,
              opt.models=opt.models)
  class(out) <- "sum.intsearch"
  
  if(verbose){
    cat("Summary interval search \n\n")
    cat(paste("show the first", first.n,"out of",nrow(out$info),"entries\n"))
    print(out$info[1:first.n,])
    cat("\n..............................")
    
    cat("\n\n Optimal parameters found are: \n\n")
    cat(paste("alpha = ",round(out$opt.alpha,digits),
              "\t",
              "lambda = ",round(out$opt.lambda,digits),
              "deviance = ",round(out$opt.error,digits)))
  }
  invisible(out)
}