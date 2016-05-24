# C-vine simulation as in Joe (2011, Chapter 7, Vine Copula Handbook),
cvinesim=function(N,param,qcondcop12,qcondcop13,qcondcop23,
                  tau2par12,tau2par13,tau2par23)
{ tau12=param[1]
  tau13=param[2]
  tau23=param[3]
  th12=tau2par12(tau12)
  th13=tau2par13(tau13)
  th23=tau2par23(tau23)
  d=3
  p=runif(N*d)
  p=matrix(p,nrow=N)
  u=matrix(0,N,d)
  u[,1]=p[,1] 
  u[,2]=qcondcop12(p[,2],p[,1], th12) 
  ttem=p[,3] 
  ttem=qcondcop23(ttem, p[,2], th23)
  ttem=qcondcop13(ttem, p[,1], th13)
  u[,3]=ttem
  u
}


rgammaShifted=function(n,shape,scale,thres) 
{ rgamma(n, shape, 1/scale) + thres }

rVineCopulaREMADA.beta=function(N,p,g,taus,qcondcop12,qcondcop13,
                           qcondcop23,tau2par12,tau2par13,tau2par23)
{ a=p/g-p
  b=(1-p)*(1-g)/g
  dat=cvinesim(N,taus,qcondcop12,qcondcop13,qcondcop23,
             tau2par12,tau2par13,tau2par23)
  u1=dat[,1]
  u2=dat[,2]
  u3=dat[,3]
  x1=qbeta(u1,a[1],b[1])
  x2=qbeta(u2,a[2],b[2])
  x3=qbeta(u3,a[3],b[3])
  n=round(rgammaShifted(N,shape=1.2,scale=100,thres=30))
  n1=round(n*x3)
  n2=n-n1
  TP=round(n1*x1)
  TN=round(n2*x2)
  FN=n1-TP
  FP=n2-TN
  list("TP"=TP,"TN"=TN,"FN"=FN,"FP"=FP)
}

rVineCopulaREMADA.norm=function(N,p,si,taus,qcondcop12,qcondcop13,
                                qcondcop23,tau2par12,tau2par13,tau2par23)
{ dat=cvinesim(N,taus,qcondcop12,qcondcop13,qcondcop23,
               tau2par12,tau2par13,tau2par23)
  u1=dat[,1]
  u2=dat[,2]
  u3=dat[,3]
  mu=log(p/(1-p))
  x1=qnorm(u1,mu[1],si[1])
  x2=qnorm(u2,mu[2],si[2])
  x3=qnorm(u3,mu[3],si[3])
  t1=exp(x1)
  t2=exp(x2)
  t3=exp(x3)
  x1=t1/(1+t1)
  x2=t2/(1+t2)
  x3=t3/(1+t3)
  n=round(rgammaShifted(N,shape=1.2,scale=100,thres=30))
  n1=round(n*x3)
  n2=n-n1
  TP=round(n1*x1)
  TN=round(n2*x2)
  FN=n1-TP
  FP=n2-TN
  list("TP"=TP,"TN"=TN,"FN"=FN,"FP"=FP)
}

