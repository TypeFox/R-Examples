jpmf.norm<-function(param,TP,FN,FP,TN,gl,mgrid,qcond,tau2par)
{ p=param[1:2]
  si=param[3:4]
  tau=param[5]
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
  prob
}


jpmf.beta<-function(param,TP,FN,FP,TN,gl,mgrid,qcond,tau2par)
{ p=param[1:2]
  g=param[3:4]
  tau=param[5]
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
  prob
}

vuong.norm=function(qcond,tau2par,param1,param2,TP,FN,FP,TN,gl,mgrid)
{ prob1=jpmf.norm(param1,TP,FN,FP,TN,gl,mgrid,qcondbvn,tau2par.bvn)
  n=length(prob1)
  prob2=jpmf.norm(param2,TP,FN,FP,TN,gl,mgrid,qcond,tau2par)
  m=log(prob2/prob1)
  z=sqrt(n)*mean(m)/sd(m)
  pvalue<-2*pnorm(-abs(z))
  result<-data.frame(round(z,digits=3),
    round(pvalue,digits=3))
  names(result)<-c("z","p.value")
  return(result)
}

vuong.beta=function(qcond,tau2par,param1,param2,TP,FN,FP,TN,gl,mgrid)
{ prob1=jpmf.norm(param1,TP,FN,FP,TN,gl,mgrid,qcondbvn,tau2par.bvn)
  n=length(prob1)
  prob2=jpmf.beta(param2,TP,FN,FP,TN,gl,mgrid,qcond,tau2par)
  m=log(prob2/prob1)
  z=sqrt(n)*mean(m)/sd(m)
  pvalue<-2*pnorm(-abs(z))
  result<-data.frame(round(z,digits=3),
                     round(pvalue,digits=3))
  names(result)<-c("z","p.value")
  return(result)
}









countermon.jpmf.beta<-function(param,TP,FN,FP,TN,gl,mgrid)
{ p=param[1:2]
  g=param[3:4]
  tau=-0.95
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
  prob
}


countermon.jpmf.norm<-function(param,TP,FN,FP,TN,gl,mgrid)
{ p=param[1:2]
  si=param[3:4]
  tau=-0.95
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
 prob
}

countermonotonicity.vuong=function(param1,param2,TP,FN,FP,TN,gl,mgrid)
{ prob1=countermon.jpmf.norm(param1,TP,FN,FP,TN,gl,mgrid)
  n=length(prob1)
  prob2=countermon.jpmf.beta(param2,TP,FN,FP,TN,gl,mgrid)
  m=log(prob2/prob1)
  z=(sqrt(n)*mean(m))/sd(m)
  pvalue<-2*pnorm(-abs(z))
  result<-data.frame(round(z,digits=3),
    round(pvalue,digits=3))
  names(result)<-c("z","p.value")
  return(result)
}


