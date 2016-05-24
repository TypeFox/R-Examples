tribinomprod=function(x1,x2,x3,TP,FN,FP,TN,perm)
{ n1=TP+FN
  n2=TN+FP
  n3=n1+n2
  if(perm==1) # 12, 13, 23|1
  { f1=dbinom(TP,size=n1,prob=x1)
    f2=dbinom(TN,size=n2,prob=x2)
    f3=dbinom((TP+FN),size=n3,prob=x3)
  } else {
    if(perm==2) # 21, 23, 13|2
    { f1=dbinom(TN,size=n2,prob=x1)
      f2=dbinom(TP+FN,size=n3,prob=x2)
      f3=dbinom(TP,size=n1,prob=x3)
      
    } else { # 31, 32, 12|3
      f1=dbinom(TP+FN,size=n3,prob=x1)
      f2=dbinom(TP,size=n1,prob=x2)
      f3=dbinom(TN,size=n2,prob=x3)
    }}
  f1*f2*f3    
}



###################################################
tvineloglik.norm<-function(param,TP,FN,FP,TN,perm,gl,mgrid,
                           qcondcop12,qcondcop13,
                          tau2par12,tau2par13)
{ p=param[1:3]
  si=param[4:6]
  tau12=param[7]
  tau13=param[8]
  if(p[1]<=0 | p[1]>=1) return(1.e10)
  if(p[2]<=0 | p[2]>=1) return(1.e10)
  if(p[3]<= 0 | p[3]>=1) return(1.e10)
  if(si[1]<=0 | si[2]<=0 | si[3]<=0) return(1.e10)
  if(tau12< -0.95 | tau12>=0.95) return(1.e10)
  if(tau13< -0.95 | tau13>=0.95) return(1.e10)
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
  -sum(log(prob))
}




tVineCopulaREMADA.norm=function(TP,FN,FP,TN,perm,gl,mgrid,
                                qcondcop12,qcondcop13,
                               tau2par12,tau2par13)
{ rTP=TP + 0.5*(TP==0)
  rFN=FN + 0.5*(FN==0)
  rFP=FP + 0.5*(FP==0)
  rTN=TN + 0.5*(TN==0)
  SE=rTP/(rTP+rFN)
  SP=rTN/(rTN+rFP)
  PR=(rTP+rFN)/(rTP+rFN+rTN+rFP)
  z=cbind(SE,SP,PR)
  if(perm==1) {sel=1:3} else {if(perm==2){sel=c(2,3,1)} else {sel=c(3,1,2)}}
  z=z[,sel]
  logitz=log(z/(1-z))
  p=apply(z,2,mean)
  si<-sqrt(apply(logitz,2,var))
  stau=cor(logitz,method="kendall")
  inipar=c(p,si,stau[1,2],stau[1,3])
  est=nlm(tvineloglik.norm,inipar,TP,FN,FP,TN,perm,gl,mgrid,
          qcondcop12,qcondcop13,
          tau2par12,tau2par13,hessian=T)
  est
}

###################################################
tvineloglik.beta<-function(param,TP,FN,FP,TN,perm,gl,mgrid,
                           qcondcop12,qcondcop13,
                          tau2par12,tau2par13)
{ p=param[1:3]
  g=param[4:6]
  tau12=param[7]
  tau13=param[8]
  if(p[1]<= 0 | p[1]>=1) return(1.e10)
  if(g[1]<= 0 | g[1]>=1) return(1.e10)
  if(p[2]<= 0 | p[2]>=1) return(1.e10)
  if(g[2]<= 0 | g[2]>=1) return(1.e10)
  if(p[3]<= 0 | p[3]>=1) return(1.e10)
  if(g[3]<= 0 | g[3]>=1) return(1.e10)
  if(tau12< -0.95 | tau12>=0.95) return(1.e10)
  if(tau13< -0.95 | tau13>=0.95) return(1.e10)
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
  -sum(log(prob))
}




tVineCopulaREMADA.beta=function(TP,FN,FP,TN,perm,gl,mgrid,
                                qcondcop12,qcondcop13,
                               tau2par12,tau2par13)
{ temp1=TP+FN
  temp2=TN+FP
  z=cbind(TP/temp1,TN/temp2,temp1/(temp1+temp2))
  if(perm==1) {sel=1:3} else {if(perm==2){sel=c(2,3,1)} else {sel=c(3,1,2)}}
  z=z[,sel]
  p<-apply(z,2,mean)
  si2<-apply(z,2,var)
  g=si2/p/(1-p)
  stau=cor(z,method="kendall")
  inipar=c(p,g,stau[1,2],stau[1,3])
  est=nlm(tvineloglik.beta,inipar,TP,FN,FP,TN,perm,gl,mgrid,
          qcondcop12,qcondcop13,
          tau2par12,tau2par13,hessian=T)
  est
}



