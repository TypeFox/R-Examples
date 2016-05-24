cv.lars <-
function(x, y, K = 10, index, 
           trace = FALSE, plot.it = TRUE, se = TRUE,type = c("lasso", "lar", "forward.stagewise", "stepwise"),
         mode=c("fraction", "step"),...)
{
  type=match.arg(type)

  if(missing(mode)){
    mode=switch(type,
      lasso="fraction",
      lar="step",
      forward.stagewise="fraction",
      stepwise="step"
      )
  }
  else  mode=match.arg(mode)
  all.folds <- cv.folds(length(y), K)
if(missing(index)){
  index=seq(from = 0, to = 1, length = 100)
  if(mode=="step"){
    fit=lars(x,y,type=type,...)
      nsteps=nrow(fit$beta)
      maxfold=max(sapply(all.folds,length))
      nsteps=min(nsteps,length(y)-maxfold)
      index=seq(nsteps)
    }
}
   residmat <- matrix(0, length(index), K)
  for(i in seq(K)) {
    omit <- all.folds[[i]]
    fit <- lars(x[ - omit,,drop=FALSE  ], y[ - omit], trace = trace, type=type,...)
    fit <- predict(fit, x[omit,  ,drop=FALSE], mode = mode, s = index
                   )$fit
    if(length(omit)==1)fit<-matrix(fit,nrow=1)
    residmat[, i] <- apply((y[omit] - fit)^2, 2, mean)
    if(trace)
      cat("\n CV Fold", i, "\n\n")
  }
  cv <- apply(residmat, 1, mean)
  cv.error <- sqrt(apply(residmat, 1, var)/K)
  object<-list(index = index, cv = cv, cv.error = cv.error,mode=mode)
  if(plot.it) plotCVLars(object,se=se)
  invisible(object)
}

