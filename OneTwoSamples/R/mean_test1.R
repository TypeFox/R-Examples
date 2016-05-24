mean_test1<-function(x, mu=0, sigma=-1, side=0){
   ## source("p_value.R")
   n<-length(x); xb<-mean(x)
   if (sigma>=0){
      z<-(xb-mu)/(sigma/sqrt(n))
      P<-p_value(pnorm, z, side=side)
      data.frame(mean=xb, df=n, Z=z, p_value=P)
   }
   else{
      t<-(xb-mu)/(sd(x)/sqrt(n))
      P<-p_value(pt, t, paramet=n-1, side=side)
      data.frame(mean=xb, df=n-1, T=t, p_value=P)
   }
}