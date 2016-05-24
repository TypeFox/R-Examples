Woolf.CI<-function(x1,n1,x2,n2,conf=0.95)
{
    if ((x1 == 0) || (x2 == 0) || (x1 == n1) || (x2== n2)) {
      return(c(0,Inf)) 
    } else {
      logtheta<-log(x1*(n2-x2)/(x2*(n1-x1)))
      CI.halflen<-qnorm(1-(1-conf)/2)*sqrt(1/x1+1/(n1-x1)+1/x2+1/(n2-x2))
      return(exp(c(logtheta-CI.halflen,logtheta+CI.halflen)))
    }
}

####################

Gart.CI<-function(x1,n1,x2,n2,conf=0.95)
{
  x1=x1+0.5
  n1=n1+1
  x2=x2+0.5
  n2=n2+1
  logtheta<-log(x1*(n2-x2)/(x2*(n1-x1)))
  CI.halflen<-qnorm(1-(1-conf)/2)*sqrt(1/x1+1/(n1-x1)+1/x2+1/(n2-x2))
  return(exp(c(logtheta-CI.halflen,logtheta+CI.halflen)))
}

####################

Agrestiind.CI<-function(x1,n1,x2,n2,conf=0.95) # Agresti independence-smoothed logit
{
    if (((x1 == 0) && (x2 == 0)) || ((x1 == n1) && (x2== n2))) {
      return(c(0,Inf))
    } else {
      x11.new=x1+2*n1*(x1+x2)/(n1+n2)^2
      x12.new=n1-x1+2*n1*(n1-x1+n2-x2)/(n1+n2)^2
      x21.new=x2+2*n2*(x1+x2)/(n1+n2)^2
      x22.new=n2-x2+2*n2*(n1-x1+n2-x2)/(n1+n2)^2
      logtheta<-log(x11.new*x22.new/(x21.new*x12.new))
      CI.halflen<-qnorm(1-(1-conf)/2)*sqrt(1/x11.new+1/x12.new+1/x21.new+1/x22.new) 
      return(exp(c(logtheta-CI.halflen,logtheta+CI.halflen)))
    }
}

#######################

Invsinh.CI<-function(x1,n1,x2,n2,phi.1,phi.2,conf=0.95) 
{
  if(phi.2==0){
    if ((x1 == 0) || (x2 == 0) || (x1 == n1) || (x2== n2)) return(c(0,Inf)) 
  } else {    
    or<-(x1+phi.1)*(n2-x2+phi.1)/(n1-x1+phi.1)/(x2+phi.1)
    SE<-2*asinh(qnorm(1-(1-conf)/2)/2*sqrt(1/(x1+phi.2)+1/(n1-x1+phi.2)+1/(x2+phi.2)+1/(n2-x2+phi.2)))
    return(exp(c(log(or)-SE,log(or)+SE))) 
  }
}

######################

MUE.func<-function(parms) # x=c(x1,n1,x2,n2)
{
  with(as.list(parms),{
    if (x1==n1)  {
      p_1u = 1
    } else {
      F_u = (x1+1)/(n1-x1)* qf(.5,2*(x1+1),2*(n1-x1) )
      p_1u = F_u/(1 + F_u)
    }
    
    if (x1==0) {
      p_1l = 0
    }  else {
      F_l= (n1-x1+1)/x1* qf(.5,2*(n1-x1+1),2*x1 )
      p_1l = 1/(1+F_l)
    }
    
    p1_mue= .5*(p_1l+p_1u)
    
    if (x2==n2) 
    {
      p_2u = 1
    } else {
      F_u = (x2+1)/(n2-x2)* qf(.5,2*(x2+1),2*(n2-x2) )
      p_2u = F_u/(1 + F_u)
    }
    
    if (x2==0) {
      p_2l = 0
    } else {
      F_l = (n2-x2+1)/x2* qf(.5,2*(n2-x2+1),2*x2 )
      p_2l = 1/(1+F_l)
    }
    
    p2_mue= .5* (p_2l+p_2u) 
    
    return(c(p1_mue*(1-p2_mue)/(p2_mue*(1-p1_mue)),p1_mue,p2_mue))
  })
}

MUE.bootstrap.prob<-function(parms,p1_mue,p2_mue) #parms=c(x1,n1,x2,n2)
{
  with(as.list(parms),{
    return(c(prob=dbinom(x1,n1,p1_mue)*dbinom(x2,n2,p2_mue),MUE.j=MUE.func(parms)[1]))
  })
}

MUE.CI<-function(x1,n1,x2,n2,conf=0.95) # MUE CI interval
{
  MUE.bootstrap<-matrix(0,(n1+1)*(n2+1),7)  
  colnames(MUE.bootstrap)<-c('x1','n1','x2','n2','prob','MUE.j','cum.prob')
  
  MUE<-MUE.func(c(x1=x1,n1=n1,x2=x2,n2=n2)) # MUE for odds ratio
  
  MUE.bootstrap[,1:4]<-as.matrix(expand.grid(x1=0:n1,n1=n1,x2=0:n2,n2=n2))
  MUE.prob<-apply(MUE.bootstrap[,1:4],1,MUE.bootstrap.prob,p1_mue=MUE[2],p2_mue=MUE[3])
  
  MUE.bootstrap[,5:6]<-t(MUE.prob)
  
  MUE.bootstrap<-MUE.bootstrap[order(MUE.bootstrap[,'MUE.j']),]
  MUE.bootstrap[,7]<-cumsum(MUE.bootstrap[,'prob'])
  
  L.i.index<-which(MUE.bootstrap[,'cum.prob']<((1-conf)/2))
  U.i.index<-which(MUE.bootstrap[,'cum.prob']<1-(1-conf)/2)    
  
  if(length(L.i.index)==0) CI.L<-0 else {
    L.i<-max(L.i.index)
    CI.L<-(MUE.bootstrap[L.i,'MUE.j']*(MUE.bootstrap[L.i+1,'cum.prob']-(1-conf)/2)+MUE.bootstrap[L.i+1,'MUE.j']*((1-conf)/2-MUE.bootstrap[L.i,'cum.prob']))/(MUE.bootstrap[L.i+1,'cum.prob']-MUE.bootstrap[L.i,'cum.prob'])
  }
  if(length(U.i.index)==((n1+1)*(n2+1)-1)) CI.U<-Inf else {
    U.i<-max(U.i.index)
    CI.U<-(MUE.bootstrap[U.i,'MUE.j']*(MUE.bootstrap[U.i+1,'cum.prob']-(1-(1-conf)/2))+MUE.bootstrap[U.i+1,'MUE.j']*((1-(1-conf)/2)-MUE.bootstrap[U.i,'cum.prob']))/(MUE.bootstrap[U.i+1,'cum.prob']-MUE.bootstrap[U.i,'cum.prob'])
  }
  return(as.vector(c(CI.L,CI.U)))
}

##############################

MOVER<-function(parms,conf=0.95)
{
  with(as.list(parms),{
    
    q1<-x1/(n1-x1)
    q2<-x2/(n2-x2)
    
    lu1<-as.vector(scoreci(x1,n1,conf.level=conf)$conf.int)
    lu2<-as.vector(scoreci(x2,n2,conf.level=conf)$conf.int)
     
    LU1<-lu1/(1-lu1)
    LU2<-lu2/(1-lu2)
    
    if(x2==n2 && x1==0) {
      
      CI.L<-0
      CI.U<-lu1[2]*(1-lu2[1])/(1-lu1[2])/lu2[1]
      
    } else if (x1==n1&&x2==0) {
      
      CI.L<-lu1[1]*(1-lu2[2])/(1-lu1[1])/lu2[2]
      CI.U<-Inf
      
    } else {
      
      CI.L<-(q1*q2-sqrt((q1*q2)^2-LU1[1]*LU2[2]*(2*q1-LU1[1])*(2*q2-LU2[2])))/(LU2[2]*(2*q2-LU2[2]))
      CI.U<-(q1*q2+sqrt((q1*q2)^2-LU1[2]*LU2[1]*(2*q1-LU1[2])*(2*q2-LU2[1])))/(LU2[1]*(2*q2-LU2[1]))
      
    }
    return(c(CI.L,CI.U))
  })
}

MOVER.CI<-function(x1,n1,x2,n2,conf=0.95) #MOVER Wilson interval
{
  
    if (((x1 == 0) && (x2 == 0)) || ((x1 == n1) && (x2== n2))) {
      
      return(c(0,Inf))
      
    } else if (x2==n2 || x1==0) {
      
      if (x2==n2 && x1==0) {
        return(MOVER(c(x1=x1,n1=n1,x2=x2,n2=n2),conf=conf))
      } else if (x2==n2 && x1!=0) {
        return(c(0,MOVER(c(x1=0,n1=n2,x2=n1-x1,n2=n1))[2]))
      } else if (x2!=n2 && x1==0){
        return(MOVER(c(x1=x1,n1=n1,x2=x2,n2=n2),conf=conf))
      }
      
    } else if (x1==n1 || x2==0) {
      
      if(x1==n1&&x2==0) {
        return(MOVER(c(x1=x1,n1=n1,x2=x2,n2=n2),conf=conf))
      } else if (x1==n1 && x2!=0) {
        return(c(1/MOVER(c(x1=0,n1=n1,x2=n2-x2,n2=n2))[2],Inf))
      } else if (x1!=n1 && x2==0) {
        return(c(MOVER(c(x1=x1,n1=n1,x2=x2,n2=n2),conf=conf)[1],Inf))
      }
      
    } else return(MOVER(c(x1=x1,n1=n1,x2=x2,n2=n2),conf=conf))

}

###########################
Cornfieldexact.CI<-function(x1,n1,x2,n2, conf = 0.95, interval = c(1e-8,1e8))
{
  
    if (((x1 == 0) && (x2 == 0)) || ((x1 == n1) && (x2== n2))) {
      LCL<-0
      UCL <- Inf
    } else if (x2==n2 || x1==0) {
      LCL<-0
      UCL<-uniroot(function(or) {sum(sapply(max(0,x1+x2-n2):x1,dFNCHypergeo, n1, n2, x1+x2, or)) - (1-conf)/2}, interval = interval)$root
    } else if (x1==n1 || x2==0) {
      LCL <- uniroot(function(or) {sum(sapply(x1:min(n1,x1+x2),dFNCHypergeo, n1, n2, x1+x2, or)) - (1-conf)/2}, interval = interval)$root
      UCL<-Inf
    } else {
      LCL <- uniroot(function(or) {sum(sapply(x1:min(n1,x1+x2),dFNCHypergeo, n1, n2, x1+x2, or)) - (1-conf)/2}, interval = interval)$root
      UCL <- uniroot(function(or) {sum(sapply(max(0,x1+x2-n2):x1,dFNCHypergeo, n1, n2, x1+x2, or)) - (1-conf)/2}, interval = interval)$root
    }
    
    return(c(LCL,UCL))
}

##########################

Cornfieldmidp.CI<-function(x1,n1,x2,n2, conf = 0.95, interval = c(1e-8,1e8))
{
    if (((x1 == 0) && (x2 == 0)) || ((x1 == n1) && (x2== n2))) {
      LCL<-0
      UCL <- Inf
    } else if (x2==n2 || x1==0) {
      LCL<-0
      UCL <- uniroot(function(or) {sum(sapply(max(0,x1+x2-n2):x1,dFNCHypergeo, n1, n2, x1+x2, or)) - dFNCHypergeo(x1,n1,n2,x1+x2,or)/2 - (1-conf)/2}, interval = interval)$root
    } else if (x1==n1 || x2==0) {
      LCL <- uniroot(function(or) {sum(sapply(x1:min(n1,x1+x2),dFNCHypergeo, n1, n2, x1+x2, or)) - dFNCHypergeo(x1,n1,n2,x1+x2,or)/2 - (1-conf)/2}, interval = interval)$root
      UCL<-Inf
    } else {
      LCL <- uniroot(function(or) {sum(sapply(x1:min(n1,x1+x2),dFNCHypergeo, n1, n2, x1+x2, or)) - dFNCHypergeo(x1,n1,n2,x1+x2,or)/2 - (1-conf)/2}, interval = interval)$root
      UCL <- uniroot(function(or) {sum(sapply(max(0,x1+x2-n2):x1,dFNCHypergeo, n1, n2, x1+x2, or)) - dFNCHypergeo(x1,n1,n2,x1+x2,or)/2 - (1-conf)/2}, interval = interval)$root
    }
    return(c(LCL,UCL))
}

###########################

or.exact.func<-function(or,x1,n1,x2,n2,conf)
{
  res<-sapply(max(0,x1+x2-n2):min(n1,x1+x2),dFNCHypergeo, n1, n2, x1+x2, or)
  index<-which(res<=dFNCHypergeo(x1,n1,n2,x1+x2,or))
  return(sum(res[index])-(1-conf))  
}

BPexact.CI<-function(x1,n1,x2,n2, conf = 0.95, interval = c(1e-8,1e8)) # Baptista-Pike exact interval
{
    if (((x1 == 0) && (x2 == 0)) || ((x1 == n1) && (x2== n2))) {
      LCL<-0
      UCL <- Inf
    } else if (x2==n2 || x1==0) {
      LCL<-0
      UCL <- uniroot(function(or) {sapply(or,or.exact.func,x1,n1,x2,n2,conf)},interval = interval)$root
    } else if (x1==n1 || x2==0) {
      LCL <- uniroot(function(or) {sapply(or,or.exact.func,x1,n1,x2,n2,conf)}, interval = interval)$root
      UCL<-Inf
    } else {     
      LCL<-uniroot(function(or) {sapply(or,or.exact.func,x1,n1,x2,n2,conf)},interval = c(interval[1], x1*(n2-x2)/x2/(n1-x1)),tol=1e-8)$root
      UCL<-uniroot(function(or) {sapply(or,or.exact.func,x1,n1,x2,n2,conf)},interval = c(x1*(n2-x2)/x2/(n1-x1) ,interval[2]),tol=1e-8)$root    
    }
    
    return(c(LCL,UCL))
}

#################################

or.midp.func<-function(or,x1,n1,x2,n2,conf)
{
  res<-sapply(max(0,x1+x2-n2):min(n1,x1+x2),dFNCHypergeo, n1, n2, x1+x2, or)
  dd<-dFNCHypergeo(x1,n1,n2,x1+x2,or)
  index<-which(res<=dd)
  return(sum(res[index])-dd/2-(1-conf))  
}

BPmidp.CI<-function(x1,n1,x2,n2, conf = 0.95, interval = c(1e-8,1e8))
{
    if (((x1 == 0) && (x2 == 0)) || ((x1 == n1) && (x2== n2))) {
      LCL<-0
      UCL <- Inf
    } else if (x2==n2 || x1==0) {
      LCL<-0
      UCL <- uniroot(function(or) {sapply(or,or.midp.func,x1,n1,x2,n2,conf)},interval = interval)$root
    } else if (x1==n1 || x2==0) {
      LCL <- uniroot(function(or) {sapply(or,or.midp.func,x1,n1,x2,n2,conf)}, interval = interval)$root
      UCL<-Inf
    } else {
      LCL<-uniroot(function(or) {sapply(or,or.midp.func,x1,n1,x2,n2,conf)},interval = c(interval[1], x1*(n2-x2)/x2/(n1-x1)),tol=1e-8)$root
      UCL<-uniroot(function(or) {sapply(or,or.midp.func,x1,n1,x2,n2,conf)},interval = c(x1*(n2-x2)/x2/(n1-x1) ,interval[2]),tol=1e-8)$root     
    }
    return(c(LCL,UCL))
}
