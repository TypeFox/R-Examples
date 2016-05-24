binomprod=function(x1,x2,TP,FN,FP,TN)
{ n1=TP+FN
  n2=TN+FP
  f1=dbinom(TP,size=n1,prob=x1)
  f2=dbinom(TN,size=n2,prob=x2)
  f1*f2
}

loglik.beta<-function(param,TP,FN,FP,TN,gl,mgrid,qcond,tau2par)
{ p=param[1:2]
  g=param[3:4]
  tau=param[5]
  if(p[1]<= 0 | p[1]>=1) return(1.e10)
  if(g[1]<= 0 | g[1]>=1) return(1.e10)
  if(p[2]<= 0 | p[2]>=1) return(1.e10)
  if(g[2]<= 0 | g[2]>=1) return(1.e10)
  if(tau< -0.95 | tau>=0.95) return(1.e10)
  u1=mgrid$x
  th=tau2par(tau)
  u2=qcond(mgrid$y,mgrid$x,th)
  a=p/g-p
  b=(1-p)*(1-g)/g
  x1=qbeta(u1,a[1],b[1])
  x2=qbeta(u2,a[2],b[2]) 
  N=length(TP)
  prob<-rep(NA,N)
  for(i in 1:N)
  { temp=binomprod(x1,x2,TP[i],FN[i],FP[i],TN[i])
    prob[i]=gl$w %*% temp %*% as.matrix(gl$w)
  }
  -sum(log(prob))
}




CopulaREMADA.beta=function(TP,FN,FP,TN,gl,mgrid,qcond,tau2par)
{ z=cbind(TP/(TP+FN),TN/(TN+FP))
  p<-apply(z,2,mean)
  si2<-apply(z,2,var)
  g=si2/p/(1-p)
  stau=cor(z,method="kendall")[1,2]
  inipar=c(p,g,stau)
  est=nlm(loglik.beta,inipar,TP,FN,FP,TN,gl,mgrid,
          qcond,tau2par,hessian=T)
  est
}


countermon.loglik.beta<-function(param,TP,FN,FP,TN,gl,mgrid)
{ p=param[1:2]
  g=param[3:4]
  tau=-0.95
  if(p[1]<= 0 | p[1]>=1) return(1.e10)
  if(g[1]<= 0 | g[1]>=1) return(1.e10)
  if(p[2]<= 0 | p[2]>=1) return(1.e10)
  if(g[2]<= 0 | g[2]>=1) return(1.e10)
  u1=mgrid$x
  th=tau2par.bvn(tau)
  u2=qcondbvn(mgrid$y,mgrid$x,th)
  a=p/g-p
  b=(1-p)*(1-g)/g
  x1=qbeta(u1,a[1],b[1])
  x2=qbeta(u2,a[2],b[2]) 
  N=length(TP)
  prob<-rep(NA,N)
  for(i in 1:N)
  { temp=binomprod(x1,x2,TP[i],FN[i],FP[i],TN[i])
    prob[i]=gl$w %*% temp %*% as.matrix(gl$w)
  }
  -sum(log(prob))
}




countermonotonicCopulaREMADA.beta=function(TP,FN,FP,TN,gl,mgrid)
{ z=cbind(TP/(TP+FN),TN/(TN+FP))
  p<-apply(z,2,mean)
  si2<-apply(z,2,var)
  g=si2/p/(1-p)
  inipar=c(p,g)
  est=nlm(countermon.loglik.beta,inipar,TP,FN,FP,TN,gl,mgrid,hessian=T)
  est
}



loglik.norm<-function(param,TP,FN,FP,TN,gl,mgrid,qcond,tau2par)
{ p=param[1:2]
  si=param[3:4]
  tau=param[5]
  if(p[1]<=0 | p[1]>=1) return(1.e10)
  if(p[2]<=0 | p[2]>=1) return(1.e10)
  if(si[1]<=0 | si[2]<=0) return(1.e10)
  if(tau< -0.95 | tau>=0.95) return(1.e10)
  mu=log(p/(1-p))
  u1=mgrid$x
  th=tau2par(tau)
  u2=qcond(mgrid$y,mgrid$x,th)
  x1=qnorm(u1,mu[1],si[1])
  x2=qnorm(u2,mu[2],si[2])
  t1=exp(x1)
  t2=exp(x2)
  x1=t1/(1+t1)
  x2=t2/(1+t2)
  N=length(TP)
  prob<-rep(NA,N)
  for(i in 1:N)
  { temp=binomprod(x1,x2,TP[i],FN[i],FP[i],TN[i])
    prob[i]=gl$w %*% temp %*% as.matrix(gl$w)
  }
  -sum(log(prob))
}




CopulaREMADA.norm=function(TP,FN,FP,TN,gl,mgrid,qcond,tau2par)
{ rTP=TP + 0.5*(TP==0)
  rFN=FN + 0.5*(FN==0)
  rFP=FP + 0.5*(FP==0)
  rTN=TN + 0.5*(TN==0)
  SE=rTP/(rTP+rFN)
  SP=rTN/(rTN+rFP)
  z=cbind(SE,SP)
  logitz=log(z/(1-z))
  p=apply(z,2,mean)
  si<-sqrt(apply(logitz,2,var))
  stau=cor(logitz,method="kendall")[1,2]
  inipar=c(p,si,stau)
  est=nlm(loglik.norm,inipar,TP,FN,FP,TN,gl,mgrid,
          qcond,tau2par,hessian=T)
  est
}

countermon.loglik.norm<-function(param,TP,FN,FP,TN,gl,mgrid)
{ p=param[1:2]
  si=param[3:4]
  tau=-0.95
  if(p[1]<=0 | p[1]>=1) return(1.e10)
  if(p[2]<=0 | p[2]>=1) return(1.e10)
  if(si[1]<=0 | si[2]<=0) return(1.e10)
  mu=log(p/(1-p))
  u1=mgrid$x
  th=tau2par.bvn(tau)
  u2=qcondbvn(mgrid$y,mgrid$x,th)
  x1=qnorm(u1,mu[1],si[1])
  x2=qnorm(u2,mu[2],si[2])
  t1=exp(x1)
  t2=exp(x2)
  x1=t1/(1+t1)
  x2=t2/(1+t2)
  N=length(TP)
  prob<-rep(NA,N)
  for(i in 1:N)
  { temp=binomprod(x1,x2,TP[i],FN[i],FP[i],TN[i])
    prob[i]=gl$w %*% temp %*% as.matrix(gl$w)
  }
  -sum(log(prob))
}

countermonotonicCopulaREMADA.norm=function(TP,FN,FP,TN,gl,mgrid)
{ rTP=TP + 0.5*(TP==0)
  rFN=FN + 0.5*(FN==0)
  rFP=FP + 0.5*(FP==0)
  rTN=TN + 0.5*(TN==0)
  SE=rTP/(rTP+rFN)
  SP=rTN/(rTN+rFP)
  z=cbind(SE,SP)
  logitz=log(z/(1-z))
  p=apply(z,2,mean)
  si<-sqrt(apply(logitz,2,var))
  inipar=c(p,si)
  est=nlm(countermon.loglik.norm,inipar,TP,FN,FP,TN,gl,mgrid,hessian=T)
  est
}

