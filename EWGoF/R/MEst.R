
#Function that computes the moment estimators of the two parameters of Weibull 
MEst<-function(x){
  if(sum(x<0)){stop(paste("Data x is not a positive sample"))}
  n = length(x)
  #Apply the log. transformation => If: X -> Weibull, then -log(X) -> EV
  y = sort(log(x))
  S=sqrt(sum((y-mean(y))^2)/(n-1))
  t = (sqrt(6)/pi)*S
  ksi = mean(y)-digamma(1)*t    #digamma(1)=-cst d'euler
  #Compute the pseudo-observation y_1, .., y_n
  y = (y-ksi)/t
  return(MEst <-list(eta=exp(ksi),beta=1/t,y=y)) 
}
