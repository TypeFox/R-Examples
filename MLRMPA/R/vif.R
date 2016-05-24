vif <-
function(y,x,lm){
  yx<-data.frame(y,x)
  colnames(yx)[1]<-c("expr")
  j<-dim(x)[2]
  r2<-numeric(j)
  vif.value<-numeric(j)
  for (i in 1:j){
    lm.new<-lm(expr~.,data=yx[,-(i+1)])
    r2[i]<-summary(lm.new)$r.squared
    vif.value[i]<-1/(1-r2[i])
  }
  VIF.value<-data.frame(VIF.value=vif.value)
  rownames(VIF.value)<-colnames(lm$model[,-1])
  return(VIF.value)
}
