mkCrisp<-function(x) {
  x<-as.factor(x)
  return(ifelse(outer(x, levels(x),"=="),1,0))
}