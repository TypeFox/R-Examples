vinepmf.norm<-function(param,TP,FN,FP,TN,perm,gl,mgrid,
                          qcondcop12,qcondcop13,qcondcop23,
                          tau2par12,tau2par13,tau2par23)
{ p=param[1:3]
  si=param[4:6]
  tau12=param[7]
  tau13=param[8]
  tau23=param[9]
  th12=tau2par12(tau12)
  th13=tau2par13(tau13)
  th23=tau2par23(tau23)
  w1=mgrid$x
  w2=mgrid$y
  w3=mgrid$z
  u1=w1
  u2=qcondcop12(mgrid$y,mgrid$x,th12)
  t=qcondcop23(w3,w2,th23)
  u3=qcondcop13(t,w1,th13)
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
  N=length(TP)
  prob<-rep(NA,N)
  for(i in 1:N)
  { temp=tribinomprod(x1,x2,x3,TP[i],FN[i],FP[i],TN[i],perm)
    prob[i]= tensor(tensor(temp,gl$w,3,1),gl$w,2,1)%*%gl$w
  }
  prob
}

vinepmf.beta<-function(param,TP,FN,FP,TN,perm,gl,mgrid,
                          qcondcop12,qcondcop13,qcondcop23,
                          tau2par12,tau2par13,tau2par23)
{ p=param[1:3]
  g=param[4:6]
  tau12=param[7]
  tau13=param[8]
  tau23=param[9]
  th12=tau2par12(tau12)
  th13=tau2par13(tau13)
  th23=tau2par23(tau23)
  w1=mgrid$x
  w2=mgrid$y
  w3=mgrid$z
  u1=w1
  u2=qcondcop12(mgrid$y,mgrid$x,th12)
  t=qcondcop23(w3,w2,th23)
  u3=qcondcop13(t,w1,th13)
  a=p/g-p
  b=(1-p)*(1-g)/g
  x1=qbeta(u1,a[1],b[1])
  x2=qbeta(u2,a[2],b[2]) 
  x3=qbeta(u3,a[3],b[3]) 
  N=length(TP)
  prob<-rep(NA,N)
  for(i in 1:N)
  { temp=tribinomprod(x1,x2,x3,TP[i],FN[i],FP[i],TN[i],perm)
    prob[i]= tensor(tensor(temp,gl$w,3,1),gl$w,2,1)%*%gl$w
  }
  prob
}

tvinepmf.norm<-function(param,TP,FN,FP,TN,perm,gl,mgrid,
                           qcondcop12,qcondcop13,
                           tau2par12,tau2par13)
{ p=param[1:3]
  si=param[4:6]
  tau12=param[7]
  tau13=param[8]
  th12=tau2par12(tau12)
  th13=tau2par13(tau13)
  w1=mgrid$x
  w2=mgrid$y
  w3=mgrid$z
  u1=w1
  u2=qcondcop12(mgrid$y,mgrid$x,th12)
  t=qcondbvn(w3,w2,0)
  u3=qcondcop13(t,w1,th13)
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
  N=length(TP)
  prob<-rep(NA,N)
  for(i in 1:N)
  { temp=tribinomprod(x1,x2,x3,TP[i],FN[i],FP[i],TN[i],perm)
    prob[i]= tensor(tensor(temp,gl$w,3,1),gl$w,2,1)%*%gl$w
  }
  prob
}

tvinepmf.beta<-function(param,TP,FN,FP,TN,perm,gl,mgrid,
                           qcondcop12,qcondcop13,
                           tau2par12,tau2par13)
{ p=param[1:3]
  g=param[4:6]
  tau12=param[7]
  tau13=param[8]
  th12=tau2par12(tau12)
  th13=tau2par13(tau13)
  w1=mgrid$x
  w2=mgrid$y
  w3=mgrid$z
  u1=w1
  u2=qcondcop12(mgrid$y,mgrid$x,th12)
  t=qcondbvn(w3,w2,0)
  u3=qcondcop13(t,w1,th13)
  a=p/g-p
  b=(1-p)*(1-g)/g
  x1=qbeta(u1,a[1],b[1])
  x2=qbeta(u2,a[2],b[2]) 
  x3=qbeta(u3,a[3],b[3]) 
  N=length(TP)
  prob<-rep(NA,N)
  for(i in 1:N)
  { temp=tribinomprod(x1,x2,x3,TP[i],FN[i],FP[i],TN[i],perm)
    prob[i]= tensor(tensor(temp,gl$w,3,1),gl$w,2,1)%*%gl$w
  }
  prob
}

vine.vuong.norm=function(qcondcop12,qcondcop13,qcondcop23,tau2par12,tau2par13,tau2par23,
                    param1,param2,TP,FN,FP,TN,perm,gl,mgrid)
{ prob1=vinepmf.norm(param1,TP,FN,FP,TN,perm,gl,mgrid,qcondbvn,qcondbvn,qcondbvn,
                  tau2par.bvn,tau2par.bvn,tau2par.bvn)
  n=length(prob1)
  prob2=vinepmf.norm(param2,TP,FN,FP,TN,perm,gl,mgrid,qcondcop12,qcondcop13,qcondcop23,
                  tau2par12,tau2par13,tau2par23)
  m=log(prob2/prob1)
  z=sqrt(n)*mean(m)/sd(m)
  pvalue<-2*pnorm(-abs(z))
  result<-data.frame(round(z,digits=3),
  round(pvalue,digits=3))
  names(result)<-c("z","p.value")
  return(result)
}

vine.vuong.beta=function(qcondcop12,qcondcop13,qcondcop23,tau2par12,tau2par13,tau2par23,
                    param1,param2,TP,FN,FP,TN,perm,gl,mgrid)
{ prob1=vinepmf.norm(param1,TP,FN,FP,TN,perm,gl,mgrid,qcondbvn,qcondbvn,qcondbvn,
                  tau2par.bvn,tau2par.bvn,tau2par.bvn)
  n=length(prob1)
  prob2=vinepmf.beta(param2,TP,FN,FP,TN,perm,gl,mgrid,qcondcop12,qcondcop13,qcondcop23,
                  tau2par12,tau2par13,tau2par23)
  m=log(prob2/prob1)
  z=sqrt(n)*mean(m)/sd(m)
  pvalue<-2*pnorm(-abs(z))
  result<-data.frame(round(z,digits=3),
  round(pvalue,digits=3))
  names(result)<-c("z","p.value")
  return(result)
}



tvine.vuong.norm=function(qcondcop12,qcondcop13,tau2par12,tau2par13,
                         param1,param2,TP,FN,FP,TN,perm,gl,mgrid)
{ prob1=vinepmf.norm(param1,TP,FN,FP,TN,perm,gl,mgrid,qcondbvn,qcondbvn,qcondbvn,
                  tau2par.bvn,tau2par.bvn,tau2par.bvn)
  n=length(prob1)
  prob2=tvinepmf.norm(param2,TP,FN,FP,TN,perm,gl,mgrid,qcondcop12,qcondcop13,
                       tau2par12,tau2par13)
  m=log(prob2/prob1)
  z=sqrt(n)*(mean(m)+1/n)/sd(m)
  pvalue<-2*pnorm(-abs(z))
  result<-data.frame(round(z,digits=3),
  round(pvalue,digits=3))
  names(result)<-c("z","p.value")
  return(result)
}

tvine.vuong.beta=function(qcondcop12,qcondcop13,tau2par12,tau2par13,
                         param1,param2,TP,FN,FP,TN,perm,gl,mgrid)
{ prob1=vinepmf.norm(param1,TP,FN,FP,TN,perm,gl,mgrid,qcondbvn,qcondbvn,qcondbvn,
                  tau2par.bvn,tau2par.bvn,tau2par.bvn)
  n=length(prob1)
  prob2=tvinepmf.beta(param2,TP,FN,FP,TN,perm,gl,mgrid,qcondcop12,qcondcop13,
                       tau2par12,tau2par13)
  m=log(prob2/prob1)
  z=sqrt(n)*(mean(m)+1/n)/sd(m)
  pvalue<-2*pnorm(-abs(z))
  result<-data.frame(round(z,digits=3),
  round(pvalue,digits=3))
  names(result)<-c("z","p.value")
  return(result)
}

tvine2.vuong.norm=function(qcondcop12,qcondcop13,tau2par12,tau2par13,
                          param1,param2,TP,FN,FP,TN,perm,gl,mgrid)
{ prob1=tvinepmf.norm(param1,TP,FN,FP,TN,perm,gl,mgrid,qcondbvn,qcondbvn,
                       tau2par.bvn,tau2par.bvn)
  n=length(prob1)
  prob2=tvinepmf.norm(param2,TP,FN,FP,TN,perm,gl,mgrid,qcondcop12,qcondcop13,
                       tau2par12,tau2par13)
  m=log(prob2/prob1)
  z=sqrt(n)*(mean(m))/sd(m)
  pvalue<-2*pnorm(-abs(z))
  result<-data.frame(round(z,digits=3),
  round(pvalue,digits=3))
  names(result)<-c("z","p.value")
  return(result)
}

tvine2.vuong.beta=function(qcondcop12,qcondcop13,tau2par12,tau2par13,
                          param1,param2,TP,FN,FP,TN,perm,gl,mgrid)
{ prob1=tvinepmf.norm(param1,TP,FN,FP,TN,perm,gl,mgrid,qcondbvn,qcondbvn,
                       tau2par.bvn,tau2par.bvn)
  n=length(prob1)
  prob2=tvinepmf.beta(param2,TP,FN,FP,TN,perm,gl,mgrid,qcondcop12,qcondcop13,
                       tau2par12,tau2par13)
  m=log(prob2/prob1)
  z=sqrt(n)*(mean(m))/sd(m)
  pvalue<-2*pnorm(-abs(z))
  result<-data.frame(round(z,digits=3),
  round(pvalue,digits=3))
  names(result)<-c("z","p.value")
  return(result)
}



















