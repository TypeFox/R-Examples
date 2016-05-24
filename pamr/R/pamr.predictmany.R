pamr.predictmany <- function(fit, newx, threshold=fit$threshold,
                             prior = fit$prior,  threshold.scale = fit$threshold.scale,
                             ...) {
  prob <-array(NA,c(length(prior),ncol(newx),length(threshold)))
  predclass <-matrix(NA,nrow=ncol(newx),ncol=length(threshold))
  
  for(i in 1:length(threshold)){
    prob[,,i] <-pamr.predict(fit,newx,threshold=threshold[i],type="posterior",...)
    predclass[,i] <-pamr.predict(fit,newx,threshold=threshold[i],type="class",...)
  }
  
  predclass <-matrix(levels(fit$y)[predclass],ncol=length((threshold)))

  return(list(prob=prob,predclass=predclass))
}













