#function that uses the fist itetation of method of moments method of moments to get the inits

giveFunction2Minimize<-function(mu,g) {
  out = function(b) (mu - (log(g*b)*(1 - b))/( log(b)*(1 - g*b)) )^2
  return(out)
}

#this function returns the suqared

giveFunction2Integrate<-function(b,g) {
  out = function(x) x^2*dmbbefd(x,b=b,g=g)
  return(out)
}

giveInits<-function(x) {
  m0<-mean(x)
  m2<-mean(x^2)
  
  #p<=1/g
  
  p0=m2 #m2 upper limit of p0
  g=1/p0
  
  #equate 1rst moment to get the mean
  myMin<-giveFunction2Minimize(mu=m0,g=g)
  b<-nlm(f=myMin,p=.1)$estimate
  
  #return a
  a=(g-1)*b/(1-g*b)
  out<-list(a=a, b=b)
  return(out)
}
