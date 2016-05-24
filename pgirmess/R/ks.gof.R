ks.gof<-function(var){
  return(ks.test(var,"pnorm",mean(var),sd(var)))
}
