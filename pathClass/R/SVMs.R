svm.fit = function(x, y, Cs, scale, DEBUG=FALSE){

  if(missing(x))     stop('No epxression matrix provided.')
  if(missing(y))     stop('No class-lables provided.')
  if(missing(Cs))    stop('No tuning parameter \'C\' provided.')
  if(missing(scale)) stop('Parameter \'scale\' must be in \'scale\', \'center\' or NULL.')
  if(length(levels(factor(y))) != 2) warning('y must have 2 levels.')
    
  scale.mean <- scale.std <- NULL

  if(!is.null(scale)){
    scale <- tolower(scale)

    if("center" %in% scale){
      x = scale(x,center=T)
      ## save centering coefficient
      scale.mean = attr(x,"scaled:center")
      names(scale.mean) = colnames(x)
    }
    
    if("scale" %in% scale){
      x = scale(x,center=F,scale=T)
      ## save scaling coefficient    
      scale.std = attr(x,"scaled:scale")
      names(scale.std) = colnames(x)
    }
  }

  ## this happens sometimes whenn all
  ## probe sets have exactely the same value
  ## and after centering everything is zero then.
  ## Due to this the scaling fails and produces
  ## NaNs:
  ## setting NaN columns to 0
  nan.cols     <- apply(x,2,function(y) all(is.nan(y)))
  x[,nan.cols] <- 0
  
  K <- kernelMatrix(vanilladot(), x)
  best.bound <- Inf

  for(C in Cs){
    if(DEBUG) cat('Trying C=',C,'\n')
    K2 = as.kernelMatrix(K + 1/C*diag(NROW(x)))
    fit.tmp = ksvm(K2, y, C=Inf, type="C-svc", shrinking=FALSE, tol=0.01,scaled=F)
    bound = spanbound(fit.tmp, K2, sign(as.numeric(y) - 1.5))
    if(bound < best.bound){
      model = fit.tmp
      best.bound = bound
      Cbest = C
    }
  }
  if(DEBUG) cat('Best C=',Cbest,'\n')
  
  ## alphaindex: The index of the resulting support vectors in the data matrix
  svs <- unlist(alphaindex(model))

  ## coef: The corresponding coefficients times the training labels.
  w <- abs(t(unlist(coef(model))) %*% x[svs,])
  
  fit <- list(fit.svm=model, w=w, K=K, C=Cbest, xsvs=x[svs,,drop=FALSE], error.bound=best.bound, scale.mean=scale.mean, scale.std=scale.std, features=colnames(x), R=NULL)	
  return(fit)
}

svm.predict = function(fit, newdata, type="response"){

  ## do the prediction only with those genes
  ## that were use for training
  newdata <- newdata[,fit$features]

  if(!is.null(fit$scale.mean))
    newdata <- scale(newdata, center=fit$scale.mean[fit$features], scale=FALSE)
  
  if(!is.null(fit$scale.std))
    newdata <- scale(newdata, center=FALSE, scale=fit$scale.std[fit$features])
  
  Ktst        <- kernelMatrix(vanilladot(), newdata, fit$xsvs[,fit$features])
  Ktst2       <- kernelMatrix(rbfdot(sigma=0.001), newdata, fit$xsvs[,fit$features])
  ident       <- which(Ktst2 == 1)		
  Ktst[ident] <- Ktst[ident] + 1/fit$C		
  alpha       <- as.matrix(unlist(coef(fit$fit.svm)))	
  yhat        <- Ktst%*%alpha - b(fit$fit.svm)

  if(type == "class") yhat <- sign(yhat)	

  return(yhat)
}
