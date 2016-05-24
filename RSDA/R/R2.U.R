R2.U <-
function(sym.var,prediction) {
  pvar<-sym.var
  pred<-prediction
  return(cor(pvar$var.data.vector[,2],pred[,2])^2)
}
