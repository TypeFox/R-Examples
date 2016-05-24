QtoJackson<-function(Q,h0, x1, x2, x3){
  Calpha<-x3*((Q*x1)^(h0)+x2)
  threshold<-pnorm(Calpha, lower.tail=FALSE)

  return(threshold)


}
