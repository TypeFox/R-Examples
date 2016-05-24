RMSE.L <-
function(sym.var,prediction) {
  pvar<-sym.var
  pred<-prediction
  nn<-pvar$N
  res<-sqrt(sum((pvar$var.data.vector[,1]-pred[,1])^2)/nn)
  #res<-sum((pvar$var.data.vector[,1]-pred[,1])^2)/nn  
  return(res)
}
