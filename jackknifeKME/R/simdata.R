simdata <-
function(n,lambda)
  {
  # n: sample size
  # lambda: parameter of Uniform dist: U(lambda, 2*lambda)
  # Returned value:
  # Y=ordered min(t,c), delta=censoring indicator
  # Pper= censoring percentage
  repeat{repeat{repeat{
  t<-rlnorm(n, meanlog=1.099, sdlog=1)
  c<-runif(n,lambda,2*lambda)
  n<-length(t)
  z<-rep(NA,n) 
  d<-rep(NA,n) 
  for(i in 1:n){
  if (t[i]<c[i]){
  z[i]<-t[i]
  d[i]<-1}
  else {z[i]<-c[i]
  d[i]<-0}
  }
  sorted<-order(z) 
  sz<-as.double(z[sorted])
  sstat<-as.integer(d[sorted])
  Pper<-sum(sstat==0)*100/length(sstat)
  if (lambda==7.53){Cper=10}
  if (lambda==4.81){Cper=20}
  if (lambda==3.48){Cper=30}
  if (lambda==2.64){Cper=40}
  if (lambda==2.04){Cper=50}
  if (lambda==1.58){Cper=60}
  if (lambda==1.20){Cper=70}
  if (lambda==0.87){Cper=80}
  if (lambda==0.55){Cper=90}
  if (sstat[n]==1) {break}
  }
  if (sstat[n-1]==0) {break}
  }
  if (Pper==Cper) {break}
  }
  list(Y=sz,delta=sstat,Pper=Pper) 
  }
