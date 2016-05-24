interval_estimate1<-function(x,sigma=-1,alpha=0.05){ 
   n<-length(x); xb<-mean(x)
   if (sigma>=0){
      tmp<-sigma/sqrt(n)*qnorm(1-alpha/2); df<-n
   }
   else{
      tmp<-sd(x)/sqrt(n)*qt(1-alpha/2,n-1); df<-n-1
   }
   data.frame(mean=xb, df=df, a=xb-tmp, b=xb+tmp)
}
