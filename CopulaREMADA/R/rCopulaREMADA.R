rgammaShifted=function (n,shape,scale,thres) 
{ rgamma(n, shape, 1/scale) + thres }

rCopulaREMADA.beta=function(N,p,g,tau,rcop,tau2par)
{ n=round(rgammaShifted(N,shape=1.2,scale=100,thres=30))
  n1=rbinom(N,size=n,prob=0.43)
  n2=n-n1
  th=tau2par(tau)
  dat=rcop(N,th)
  u1=dat[,1]
  u2=dat[,2]
  a=p/g-p
  b=(1-p)*(1-g)/g
  x1=qbeta(u1,a[1],b[1])
  x2=qbeta(u2,a[2],b[2])
  TP=round(n1*x1)
  TN=round(n2*x2)
  FN=n1-TP
  FP=n2-TN
  list("TP"=TP,"TN"=TN,"FN"=FN,"FP"=FP)
}

rCopulaREMADA.norm=function(N,p,si,tau,rcop,tau2par)
{ n=round(rgammaShifted(N,shape=1.2,scale=100,thres=30))
  n1=rbinom(N,size=n,prob=0.43)
  n2=n-n1
  th=tau2par(tau)
  dat=rcop(N,th)
  u1=dat[,1]
  u2=dat[,2]
  mu=log(p/(1-p))
  x1=qnorm(u1,mu[1],si[1])
  x2=qnorm(u2,mu[2],si[2])
  t1=exp(x1)
  t2=exp(x2)
  x1=t1/(1+t1)
  x2=t2/(1+t2)
  TP=round(n1*x1)
  TN=round(n2*x2)
  FN=n1-TP
  FP=n2-TN
  list("TP"=TP,"TN"=TN,"FN"=FN,"FP"=FP)
}