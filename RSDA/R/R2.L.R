R2.L <-
function(sym.var,prediction) {
  pvar<-sym.var
  pred<-prediction
  return(cor(pvar$var.data.vector[,1],pred[,1])^2)
}
