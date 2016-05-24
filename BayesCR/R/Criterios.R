LogVerosCens<-function(cc,y,mu,sigmae,lambda,nu,type="Normal",cens="right")
{
  if(cens=="left"){cens <- "1"}
  if(cens=="right"){cens <- "2"}
  
  m <- length(y)
  ver<- matrix(0,m,1)
  sigma2 <- sigmae
  censtype <- cens
  auxy<- matrix(0,m,1)
  
  if(type=="Normal") 
  {
    auxy <- (y-mu)/sqrt(sigmae)
    if (cens=="2")
    {
      ver[cc==0] <- dnorm(auxy[cc==0])/sqrt(sigmae)
      ver[cc==1] <- pnorm(-auxy[cc==1])
    }
    if (cens=="1")
    {
      ver[cc==0]<-dnorm(auxy[cc==0])/sqrt(sigmae)
      ver[cc==1]<-pnorm(auxy[cc==1])
    }
    
  }    
  if(type=="T") 
  {      
    auxy <- (y-mu)/sqrt(sigmae)
    if (cens=="2")
    {
      ver[cc==0]<-dt(auxy[cc==0],df=nu)/sqrt(sigmae)
      ver[cc==1]<-pt(-auxy[cc==1],df=nu)
    }
    if (cens=="1")
    {
      ver[cc==0]<-dt(auxy[cc==0],df=nu)/sqrt(sigmae)
      ver[cc==1]<-pt(auxy[cc==1],df=nu)
    }  
    
  }              
  if(type=="Slash") 
  {  
    auxy <- (y-mu)/sqrt(sigmae)
    if (cens=="2")
    {
      ver[cc==0]<-dSlash(auxy[cc==0],0,1,nu)/sqrt(sigmae)
      ver[cc==1]<-AcumSlash(-auxy[cc==1],0,1,nu)
    }
    if (cens=="1")
    {
      ver[cc==0]<-dSlash(auxy[cc==0],0,1,nu)/sqrt(sigmae)
      ver[cc==1]<-AcumSlash(auxy[cc==1],0,1,nu)
    }
   
  }  
  if(type=="NormalC") 
  {    
    
    auxy <- (y-mu)/sqrt(sigmae)
    if (cens=="2")
    {
      ver[cc==0]<-dNormalC(auxy[cc==0],0,1,nu)/sqrt(sigmae)
      ver[cc==1]<-AcumNormalC(-auxy[cc==1],0,1,nu)
    }
    if (cens=="1")
    {
      ver[cc==0]<-dNormalC(auxy[cc==0],0,1,nu)/sqrt(sigmae)
      ver[cc==1]<-AcumNormalC(auxy[cc==1],0,1,nu)
    }
    
  }
  
  if(type=="SN") 
  {
    ver[cc==0]<-pdfSNI(y[cc==0],mu[cc==0],sigma2,lambda,0,type="SN")
    if(cens==1)
    {
      ver[cc==1]<-cdfSNI(y[cc==1],mu[cc==1],sigma2,lambda,0,type="SN")
    }
    if(cens=="2")
    {  
      ver[cc==1]<-(1-cdfSNI(y[cc==1],mu[cc==1],sigma2,lambda,0,type="SN"))
    }
  }
  
  if(type=="ST") 
  {
    ver[cc==0]<-pdfSNI(y[cc==0],mu[cc==0],sigma2,lambda,nu,type="ST")
    if(cens==1)
    {
      ver[cc==1]<-cdfSNI(y[cc==1],mu[cc==1],sigma2,lambda,nu,type="ST")
    }
    if(cens=="2")
    {   
      ver[cc==1]<-1-cdfSNI(y[cc==1],mu[cc==1],sigma2,lambda,nu,type="ST")
    }  
  }
  if(type=="SSL") 
  {
    for (j in 1:m)
    {
      slash <- cdfSNI(y[j],mu[j],sigma2,lambda,nu,type="SSL")
      if(cc[j]==0){ver[j]<-slash$pdf}
      if(sum(cc)>0){
        if(cc[j]==1) 
        {
          if(cens=="1")
          {   
            ver[j]<-slash$cdf 
          }
          if(cens=="2")
          {   
            ver[j]<-1 - slash$cdf
          }  
        }
      } 
    }
  }
  return(ver)
}

criterios<-function(cc,y,espac=20,cadeia,type="T", cens="2", p=p,influence){
  m<-length(y)
  if (type == "Normal")
  {
    ver<-LogVerosCens(cc,y,apply(cadeia$mu,2,mean),mean(cadeia$sigma2),0,1,type=type, cens=cens)
    n.iter<-length(cadeia$sigma2)/espac
    iter <- n.iter
    Loglikaux<-matrix(0,m,n.iter)
    CPOaux<-matrix(0,m,n.iter)
    for(k in 1:n.iter){
      i<-espac*k
      fss<-cadeia$mu[k,]
      sigma2es<-cadeia$sigma2[k]
      Loglikaux[,k]<- LogVerosCens(cc,y,fss,sigma2es,0,50,type=type, cens=cens)
      CPOaux[,k]<-1/Loglikaux[,k]
    }
    Np<-p+1
  }
  
  if (type == "T")
  {
    ver<-LogVerosCens(cc,y,apply(cadeia$mu,2,mean),mean(cadeia$sigma2),0,mean(cadeia$nu),type=type, cens=cens)
    n.iter<-length(cadeia$sigma2)/espac
    iter <- n.iter
    Loglikaux<-matrix(0,m,n.iter)
    CPOaux<-matrix(0,m,n.iter)
    for(k in 1:n.iter){
      i<-espac*k
      fss<-cadeia$mu[k,]
      sigma2es<-cadeia$sigma2[k]
      nuss<-cadeia$nu[k]
      Loglikaux[,k]<- LogVerosCens(cc,y,fss,sigma2es,0,nuss,type=type, cens=cens)
      CPOaux[,k]<-1/Loglikaux[,k]
    }
    Np<-p+2
  }
  
  if (type == "Slash")
  {
    ver<-LogVerosCens(cc,y,apply(cadeia$mu,2,mean),mean(cadeia$sigma2),0,mean(cadeia$nu),type=type, cens=cens)
    n.iter<-length(cadeia$sigma2)/espac
    iter <- n.iter
    Loglikaux<-matrix(0,m,n.iter)
    CPOaux<-matrix(0,m,n.iter)
    for(k in 1:n.iter){
      i<-espac*k
      fss<-cadeia$mu[k,]
      sigma2es<-cadeia$sigma2[k]
      nuss<-cadeia$nu[k]
      Loglikaux[,k]<- LogVerosCens(cc,y,fss,sigma2es,0,nuss,type=type, cens=cens)
      CPOaux[,k]<-1/Loglikaux[,k]
    }
    Np<-p+2
  }
  
  if (type == "NormalC")
  {
    ver<-LogVerosCens(cc,y,apply(cadeia$mu,2,mean),mean(cadeia$sigma2),0,c(mean(cadeia$nu),mean(cadeia$rho)),type=type, cens=cens)
    n.iter<-length(cadeia$sigma2)/espac
    iter <- n.iter
    Loglikaux<-matrix(0,m,n.iter)
    CPOaux<-matrix(0,m,n.iter)
    for(k in 1:n.iter){
      i<-espac*k
      fss<-cadeia$mu[k,]
      sigma2es<-cadeia$sigma2[k]
      nuss<-c(cadeia$nu[k],cadeia$rho[k])
      Loglikaux[,k]<- LogVerosCens(cc,y,fss,sigma2es,0,nuss,type=type, cens=cens)
      CPOaux[,k]<-1/Loglikaux[,k]
    }
    Np<-p+3
  }
  
  if (type == "SN")
  {
    ver<-LogVerosCens(cc,y,apply(cadeia$mu,2,mean),mean(cadeia$sigma2),mean(cadeia$lambda),0,type="SN",cens=cens)
    n.iter<-length(cadeia$sigma2)/espac
    iter <- n.iter
    Loglikaux<-matrix(0,m,n.iter)
    CPOaux<-matrix(0,m,n.iter)
    for(k in 1:n.iter)
    {
      i <- espac*k
      fss<-cadeia$mu[k,]
      sigma2es<-cadeia$sigma2[k]
      lambdas<- cadeia$lambda[k]
      Loglikaux[,k]<- LogVerosCens(cc,y,fss,sigma2es,lambdas,50,type="SN",cens=cens)
      CPOaux[,k]<-1/Loglikaux[,k]
    }
    Np<-p+2
  }
  
  if (type == "ST")
  {
    ver<-LogVerosCens(cc, y,apply(cadeia$mu,2,mean),mean(cadeia$sigma2),mean(cadeia$lambda),mean(cadeia$nu),type="ST",cens=cens)
    n.iter<-length(cadeia$sigma2)/espac
    iter <- n.iter
    Loglikaux<-matrix(0,m,iter)
    CPOaux<-matrix(0,m,iter)
    for(k in 1:iter)
    {
      i<-espac*k
      fss<-cadeia$mu[k,]
      sigma2es<-cadeia$sigma2[k]
      lambdas<- cadeia$lambda[k]
      nuss<-cadeia$nu[k]
      Loglikaux[,k]<- LogVerosCens(cc,y,fss,sigma2es,lambdas,nuss,type="ST",cens=cens)
      CPOaux[,k]<-1/Loglikaux[,k]
    }
    Np<-p+3
  }
  if (type == "SSL")
  {
    ver<-LogVerosCens(cc,y,apply(cadeia$mu,2,mean),mean(cadeia$sigma2),mean(cadeia$lambda),mean(cadeia$nu),type="SSL",cens=cens)
    n.iter<-length(cadeia$sigma2)/espac
    iter <- n.iter
    Loglikaux<-matrix(0,m,iter)
    CPOaux<-matrix(0,m,iter)
    for(k in 1:iter)
    {
      i<-espac*k
      fss<-cadeia$mu[k,]
      sigma2es<-cadeia$sigma2[k]
      lambdas<- cadeia$lambda[k]
      nuss<-cadeia$nu[k]
      Loglikaux[,k]<- LogVerosCens(cc,y,fss,sigma2es,lambdas,nuss,type="SSL",cens=cens)
      CPOaux[,k]<-1/Loglikaux[,k]
    }
    Np<-p+3
  }
  
  
  mean.log.p <- apply(log(Loglikaux),1,mean)
  mean.p <- apply(Loglikaux,1,mean)
  var.log <- matrix(0,m,iter)
  for(k in 1:iter)
  {
    var.log[,k] <- (log(Loglikaux[,k])-mean.log.p)^2
  }
  var.log.p <- apply(var.log,1,sum)/(iter-1)
  
  CPO<-sum(log(1/(apply(CPOaux,1,mean))))
  pdic <- 2*(sum(log(ver))-sum(mean.log.p))
  DIC <- 2*pdic - 2*sum(log(ver))
  EAIC<- -2*sum(log(ver))+2*Np
  EBIC<- -2*sum(log(ver))+log(m)*Np
  WAIC1 <- -2*(sum(log(mean.p))-2*sum(log(mean.p)-mean.log.p))
  WAIC2 <- -2*(sum(log(mean.p))-sum(var.log.p))
  
  if(influence=="FALSE")
  {
    return(list(CPO=CPO,DIC=DIC,EAIC=EAIC,EBIC=EBIC,WAIC1=WAIC1,WAIC2=WAIC2))
  }  
  else{
    aa<- exp(log(apply(CPOaux,1,mean))+ log(Loglikaux))
    IL<- apply(log(aa),1,mean)
    JL<- apply((aa-1)*log(aa),1,mean)
    LL<- apply(abs(aa-1),1,mean)
    CHL<- apply((aa-1)^2,1,mean)
    
    return(list(CPO=CPO,DIC=DIC,EAIC=EAIC,EBIC=EBIC,WAIC1=WAIC1,WAIC2=WAIC2,KL=IL,JDist=JL,LDist=LL,ChiDist=CHL))
  }
}