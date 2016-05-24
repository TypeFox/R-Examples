boot.yhat<-function(data,indices,lmOut,regrout0){
  data<-data[indices,]
  blmOut<-lm(formula=lmOut$call$formula,data=data)
  regrout<-calc.yhat(blmOut)
  pm<-as.vector(regrout$PredictorMetrics[-nrow(regrout$PredictorMetrics),])
  apsm<-as.vector(as.matrix(regrout$APSRelatedMetrics[-nrow(regrout$APSRelatedMetrics),-2]))
  pdm<-as.vector(as.matrix(regrout$PairedDominanceMetrics))
  tau<-vector(length=ncol(regrout$PredictorMetrics))
  for (i in 1:length(tau)){
    s1<-(regrout0$PredictorMetrics[1:(nrow(regrout0$PredictorMetrics)-1),i])
    s2<-(regrout$PredictorMetrics[1:(nrow(regrout$PredictorMetrics)-1),i])
    tau[i]<-cor.test(s1,s2,method="kendall",exact=FALSE)$estimate
  }
  c(pm,pdm,apsm,tau)
}