mpbt <-
  function(data, method, type, k) {

 nn.cl=function(y1, s1, y2, s2){
  n=length(y1)
  s1.sq=s1^2
  s2.sq=s2^2
    
  y1[is.na(y1)==1]=0
  y2[is.na(y2)==1]=0
  s1.sq[is.na(s1.sq)==1]=10^10
  s2.sq[is.na(s2.sq)==1]=10^10
  
  dat.obs=data.frame(cbind(y1, y2, s1.sq, s2.sq))
  colnames(dat.obs)=c("y1", "y2", "s1.sq", "s2.sq")
  
  
  fit1 = mvmeta(y1,s1.sq, data=dat.obs,method="reml")
  fit2 = mvmeta(y2,s2.sq, data=dat.obs,method="reml")
  tau1.sq.est=fit1$Psi[1,1]
  tau2.sq.est=fit2$Psi[1,1]
  
  A=A.star=array(NA, c(n, 2,2))
  SND=P=SND.star=P.star=array(NA, c(n, 2))
  s1.sq=dat.obs$s1.sq
  s2.sq=dat.obs$s2.sq
  y1=dat.obs$y1
  y2=dat.obs$y2
  y=cbind(y1, y2)
  
  for (i in 1:n){
    A[i,,]=matrix(c((s1.sq[i]+tau1.sq.est)^(-0.5), 0, 0, (s2.sq[i]+tau2.sq.est)^(-0.5)), ncol=2)
    SND[i,]=y[i,]%*%A[i,,]
    P[i,]=c(1,1)%*%A[i,,]
  }
  
  SND1=SND[,1]
  SND2=SND[,2]
  P1=P[,1]
  P2=P[,2]
  
  
  est1.ftn=function(SND1, SND2, P1, P2){
    b1=(sum(SND1)-sum(SND1*P1))/(sum(P1)-sum(P1^2))
    a1=(sum(SND1)-b1*sum(P1))/n
    b2=(sum(SND2)-sum(SND2*P2))/(sum(P2)-sum(P2^2))
    a2=(sum(SND2)-b1*sum(P2))/n
    est=c(a1,b1, a2, b2)
    return(est)
  }
  
  est2.ftn=function(SND1, SND2, P1, P2){
    b1=sum(SND1*P1)/sum(P1^2)
    b2=sum(SND2*P2)/sum(P2^2)
    est=c(b1,b2)
    return(est)
  }
 
  est1=est1.ftn(SND1, SND2, P1, P2)
  est2=est2.ftn(SND1, SND2, P1, P2)
  

 
  ###### CL ##############################
  loglik.CL=function(t,SND1, SND2, P1, P2){
    n=length(SND1)
    mu1=t[1]+t[2]*P1
    mu2=t[3]+t[4]*P2
    loglik=-0.5*sum((SND1-mu1)^2)-0.5*sum((SND2-mu2)^2)
    return(loglik)
  }
  
  loglik.CL0=function(tt,SND1, SND2, P1, P2){
    t=c(0, tt[1], 0, tt[2])
    result=loglik.CL(t, SND1, SND2, P1, P2)
    return(result)
  }
  
  ##### score tests ###########################
  sc=sc0=array(NA, c(4,n))
  sc0.sq=sc.sq=hc0=hc=array(NA, c(4,4,n))
  
  for (j in 1:n){
     
    sc[,j]=as.matrix(jacobian(loglik.CL, x=est1, SND1=SND1[j], P1=P1[j], SND2=SND2[j], P2=P2[j]), ncol=4)
    sc0[,j]=as.matrix(jacobian(loglik.CL, x=c(0, est2[1], 0, est2[2]), SND1=SND1[j], P1=P1[j], SND2=SND2[j], P2=P2[j]), ncol=4)
    
  }
  
  JC0=cov(t(sc0))
  esc=apply(sc0, 1, mean)
  
  sct=n*t(esc)%*%solve(JC0)%*%esc
  
  
  ######
  u.sc.gg=n*esc[c(1,3)]
  I0.sc=-as.matrix(hessian(loglik.CL, x=c(0, est2[1], 0, est2[2]), SND1=SND1, P1=P1, SND2=SND2, P2=P2), ncol=4)
  I1.sc=cov(t(sc0))
  
  sigma.sc=solve(I0.sc)%*%I1.sc%*%solve(I0.sc)
  sigma.sc.gg=sigma.sc[c(1,3), c(1,3)]
  sigma.sc.inv=solve(sigma.sc)
  sigma.sc.inv.gg=sigma.sc.inv[c(1,3), c(1,3)]
  
  I0.sc.inv=solve(I0.sc)
  I0.sc.inv.gg=I0.sc.inv[c(1,3), c(1,3)]
  
  sc.ec=1/n*t(u.sc.gg)%*%I0.sc.inv.gg%*%solve(sigma.sc.gg)%*%I0.sc.inv.gg%*%u.sc.gg
  sc.mb=1/n*t(u.sc.gg)%*%I0.sc.inv.gg%*%u.sc.gg
  mylambda.sc=eigen(solve(I0.sc.inv.gg)%*%sigma.sc.gg)$values
  bar.lambda.sc=mean(mylambda.sc)
  sc.mb.a=sc.mb/bar.lambda.sc
  pv=pchisq(q=sc.mb.a, df=2, lower.tail=FALSE)
  return(list(mpbt.TS=sc.mb.a, pv=pv))
}

    if (missing(data)){data=NULL}
    if (is.null(data)){
       stop("The dataset must be specified.")
    }
 
    if (missing(type)){type=NULL}
    if (is.null(type)){
       stop('The type of outcome must be specified. Please input "continuous" or "binary". ')
    }    
  
    if (missing(k)){k=NULL}
    if (is.null(k)){
       stop("The number of outcomes must be specified.")
    } 
    if (k>2){
       stop("The method for MMA with more than 2 outcomes is currently under development. ")
    } 

    if (missing(method)){method=NULL}
    if (is.null(method)){
       stop("A method must be specified. ")
    }

    if (type=="binary"){
       stop("The method for meta-analysis with binary outcome is currently
under development. ")
}

    if (type=="continuous"){
      y1=data$y1
      s1=data$s1
      y2=data$y2
      s2=data$s2

      if(method=="nn.cl"){
        fit=nn.cl(y1, s1, y2, s2)
      }
    
      if ((method%in%c("nn.cl"))!=1){
        stop("The input method is not available for continuous data")
      }
    }  

    res=list(type=type, k=k, method=method, mpbt.TS=fit$mpbt.TS, mpbt.pv=fit$pv)
    class(res) = c("mpbt")
    return(res)
  }
