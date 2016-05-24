interval_estimate3<-function(x,sigma=-1,alpha=0.05){ 
   n<-length(x); xb<-mean(x)
   if (sigma>=0)
      tmp<-sigma/sqrt(n)*qnorm(1-alpha/2)
   else
      tmp<-sd(x)/sqrt(n)*qnorm(1-alpha/2)
   data.frame(mean=xb, a=xb-tmp, b=xb+tmp)
}