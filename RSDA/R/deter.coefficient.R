deter.coefficient <-
function(sym.var,prediction) {
  pvar<-sym.var
  nn<-pvar$N
  centers.pvar<-rep(0,nn)
  pred<-prediction
  centers.pred<-rep(0,nn)
  centers.pvar<-(pvar$var.data.vector[,1]+pvar$var.data.vector[,2])/2
  centers.pred<-(pred[,1]+pred[,2])/2
  coef<-sum((centers.pred-mean(centers.pvar))^2)/sum((centers.pvar-mean(centers.pvar))^2)
  return(coef)  
}
