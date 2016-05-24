
#Function that computes the least squares estimators of the two parameters of Weibull 
LSEst<-function(x){
  if(sum(x<0)){stop(paste("Data x is not a positive sample"))}
  n=length(x)
  #Apply the log. transformation => If: X -> Weibull, then -log(X) -> EV
  lv = -log(x)
  y = sort(lv)
  l=seq(1,n)
  p=(l-.5)/n
  c=-log(-log(p))
  cb=mean(c)
  yb=mean(y)
  b=sum((y-yb)*(c-cb))/sum((c-cb)^2)
  ksi=yb-b*cb
  #Compute the pseudo-observation y_1, .., y_n
  y = -(y-ksi)/b
  return(LSEst<-list(eta=exp(-ksi),beta=1/b,y=y))
}

