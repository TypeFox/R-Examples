RMSE.U <-
function(sym.var,prediction) {
  pvar<-sym.var
  pred<-prediction
  nn<-pvar$N
  res<-sqrt(sum((pvar$var.data.vector[,2]-pred[,2])^2)/nn)
  #res<-sum((pvar$var.data.vector[,2]-pred[,2])^2)/nn
  return(res)
}
