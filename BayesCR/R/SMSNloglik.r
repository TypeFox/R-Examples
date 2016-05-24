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
      if(sum(cc)>0){ver[cc==1] <- pnorm(-auxy[cc==1])}
    }
    if (cens=="1")
    {
      ver[cc==0]<-dnorm(auxy[cc==0])/sqrt(sigmae)
      if(sum(cc)>0){ver[cc==1]<-pnorm(auxy[cc==1])}
    }
    
  }    
  if(type=="T") 
  {      
    auxy <- (y-mu)/sqrt(sigmae)
    if (cens=="2")
    {
      ver[cc==0]<-dt(auxy[cc==0],df=nu)/sqrt(sigmae)
      if(sum(cc)>0){ver[cc==1]<-pt(-auxy[cc==1],df=nu)}
    }
    if (cens=="1")
    {
      ver[cc==0]<-dt(auxy[cc==0],df=nu)/sqrt(sigmae)
      if(sum(cc)>0){ver[cc==1]<-pt(auxy[cc==1],df=nu)}
    }  
    
  }              
  if(type=="Slash") 
  {  
    auxy <- (y-mu)/sqrt(sigmae)
    if (cens=="2")
    {
      ver[cc==0]<-dSlash(auxy[cc==0],0,1,nu)/sqrt(sigmae)
      if(sum(cc)>0){ver[cc==1]<-AcumSlash(-auxy[cc==1],0,1,nu)}
    }
    if (cens=="1")
    {
      ver[cc==0]<-dSlash(auxy[cc==0],0,1,nu)/sqrt(sigmae)
      if(sum(cc)>0){ver[cc==1]<-AcumSlash(auxy[cc==1],0,1,nu)}
    }
    
  }  
  if(type=="NormalC") 
  {    
    
    auxy <- (y-mu)/sqrt(sigmae)
    if (cens=="2")
    {
      ver[cc==0]<-dNormalC(auxy[cc==0],0,1,nu)/sqrt(sigmae)
      if(sum(cc)>0){ver[cc==1]<-AcumNormalC(-auxy[cc==1],0,1,nu)}
    }
    if (cens=="1")
    {
      ver[cc==0]<-dNormalC(auxy[cc==0],0,1,nu)/sqrt(sigmae)
      if(sum(cc)>0){ver[cc==1]<-AcumNormalC(auxy[cc==1],0,1,nu)}
    }
    
  }
  
  if(type=="SN") 
  {
    ver[cc==0]<-pdfSNI(y[cc==0],mu[cc==0],sigma2,lambda,0,type="SN")
    if(cens==1)
    {
      if(sum(cc)>0){ver[cc==1]<-cdfSNI(y[cc==1],mu[cc==1],sigma2,lambda,0,type="SN")}
    }
    if(cens=="2")
    {  
      if(sum(cc)>0){ver[cc==1]<-(1-cdfSNI(y[cc==1],mu[cc==1],sigma2,lambda,0,type="SN"))}
    }
  }
  
  if(type=="ST") 
  {
    ver[cc==0]<-pdfSNI(y[cc==0],mu[cc==0],sigma2,lambda,nu,type="ST")
    if(cens=="1")
    {
      if(sum(cc)>0){ver[cc==1]<-cdfSNI(y[cc==1],mu[cc==1],sigma2,lambda,nu,type="ST")}
    }
    if(cens=="2")
    {   
      if(sum(cc)>0){ver[cc==1]<-1-cdfSNI(y[cc==1],mu[cc==1],sigma2,lambda,nu,type="ST")}
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

