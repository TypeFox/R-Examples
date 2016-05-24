summary.bs.plmm <-
function(object, probs=c(0.005,0.025,0.05,0.95,0.975,0.995), ...){
  BS_mean<-colMeans(object[[1]])
  BC_sd<-apply(object[[1]], 2, sd)
  BS_q<-apply(object[[1]], 2, quantile, probs=probs)
  
  #VC_mean=mean(object[[1]][,len])
  #VC_sd=sd(object[[1]][,len])
  #VC_q=quantile(object[[1]][,len], probs=probs)
  #coef=list(mean=coef_mean, sd=coef_sd, quantiles=coef_q)
  #var.comp=list(mean=VC_mean, sd=VC_sd, quantiles=VC_q)
  res<-list(mean=BS_mean, sd=BC_sd, quantiles=BS_q)
  #res=list(mean=BS_mean, sd=BS_sd, quantiles=BS_quantiles)
  class(res)<-"summary.bs.plmm"
  return(res)
}
