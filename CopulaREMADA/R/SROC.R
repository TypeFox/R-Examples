SROC.beta=function(param,dcop,qcondcop,tau2par,TP,FN,FP,TN,points=TRUE,curves=TRUE)
{ p=param[1:2]
  g=param[3:4]
  tau=param[5]
  th=tau2par(tau)
  a=p/g-p
  b=(1-p)*(1-g)/g
  x1=seq(0.001,0.999,length=100)
  x2=seq(0.001,0.999,length=100)
  du=dbeta(x1,a[1],b[1])
  dv=dbeta(x2,a[2],b[2])
  u=pbeta(x1,a[1],b[1])
  v=pbeta(x2,a[2],b[2])
  res<-matrix(NA,length(u),length(v))
  for (i in 1:length(u))
  { for (j in 1:length(v))
  { res[i,j]<-du[i]*dv[j]*dcop(u[i],v[j],th)
  }
  }
  contour(x1,x2,res,nlevels=20,col=1,drawlabels=FALSE,
          xlab="Sensitivity",ylab="Specificity",lty=2,
          xlim=c(0,1),ylim=c(0,1))
  points(p[1],p[2],pch=15,cex = 2)
  if(points)
  { z=cbind(TP/(TP+FN),TN/(TN+FP))
    points(z)
  }
  if(curves)
  { quant=c(0.99,0.5,0.01)
    linesType=c(3,1,3)
    for(j in 1:3)
    { ustar=qcondcop(quant[j],v,th)
      x2star=qbeta(ustar,a[1],b[1])
      vstar=qcondcop(quant[j],u,th)
      x1star=qbeta(vstar,a[2],b[2])
      lines(x2,x1star,col=2,lty=linesType[j],lwd=1.5)
      lines(x2star,x1,col=3,lty=linesType[j],lwd=1.5)
    }
  }
}



SROC.norm=function(param,dcop,qcondcop,tau2par,TP,FN,FP,TN,points=TRUE,curves=TRUE)
{ p=param[1:2]
  si=param[3:4]
  tau=param[5]
  mu=log(p/(1-p))
  th=tau2par(tau)
  x1=seq(-5,5,length=100)
  x2=seq(-5,5,length=100)
  du=dnorm(x1,mu[1],si[1])
  dv=dnorm(x2,mu[2],si[2])
  u=pnorm(x1,mu[1],si[1])
  v=pnorm(x2,mu[2],si[2])
  res<-matrix(NA,length(u),length(v))
  for (i in 1:length(u))
  { for (j in 1:length(v))
  { res[i,j]<-du[i]*dv[j]*dcop(u[i],v[j],th)
  }
  }
  contour(x1,x2,res,nlevels=20,col=1,drawlabels=FALSE,lty=2,
          xlab="logit(Sensitivity)",ylab="logit(Specificity)",
          xlim=c(-5,5),ylim=c(-5,5))
  points(mu[1],mu[2],pch=15,cex = 2)
  if(points)
  { rTP=TP + 0.5*(TP==0)
    rFN=FN + 0.5*(FN==0)
    rFP=FP + 0.5*(FP==0)
    rTN=TN + 0.5*(TN==0)
    SE=rTP/(rTP+rFN)
    SP=rTN/(rTN+rFP)
    z=cbind(SE,SP)
    logitz=log(z/(1-z))
    points(logitz[,1],logitz[,2])
  }
  if(curves)
  { quant=c(0.99,0.5,0.01)
    linesType=c(3,1,3)
    for(j in 1:3)
    { ustar=qcondcop(quant[j],v,th)
      x2star=qnorm(ustar,mu[1],si[1])
      vstar=qcondcop(quant[j],u,th)
      x1star=qnorm(vstar,mu[2],si[2])
      lines(x2,x1star,col=2,lty=linesType[j],lwd=1.5)
      lines(x2star,x1,col=3,lty=linesType[j],lwd=1.5)
    }
  }  
}

SROC=function(param.beta,param.normal,TP,FN,FP,TN)
{ p.beta=param.beta[1:2]
  g=param.beta[3:4]
  p.normal=param.normal[1:2]
  si=param.normal[3:4]
  mu=log(p.normal/(1-p.normal))
  th.normal=-0.99
  a=p.beta/g-p.beta
  b=(1-p.beta)*(1-g)/g
  x=seq(0.001,0.999,length=100)
  u.beta=pbeta(x,a[1],b[1])
  vstar.beta=qcondbvn(0.5,u.beta,th.normal)
  x1star.beta=qbeta(vstar.beta,a[2],b[2])
  
  logitx=log(x/(1-x))
  u.normal=pnorm(logitx,mu[1],si[1])
  vstar.normal=qcondbvn(0.5,u.normal,th.normal)
  x1star.normal=qnorm(vstar.normal,mu[2],si[2])
  x1star.normal=exp(x1star.normal)
  x1star.normal=x1star.normal/(1+x1star.normal)
  
  z=cbind(TP/(TP+FN),TN/(TN+FP))
  
  plot(z,xlab="Sensitivity",ylab="Specificity",xlim=c(0,1),ylim=c(0,1),type="p")
  points(p.normal[1],p.normal[2],pch=15,cex = 1.5)
  lines(x,x1star.normal,col=1,lty=1)
  points(p.beta[1],p.beta[2],pch=18,cex = 2,col=2)
  lines(x,x1star.beta,col=2,lty=2)
  legend(0,0.3,c('',''),lty=1:2,col=1:2,cex=0.7,bty="n")
  legend(0.1,0.3,c('Normal','Beta'),pch=c(15,18),col=1:2,cex=0.6,bty="n")
  
}




