interval_estimate4<-function(x, sigma=-1, side=0, alpha=0.05){ 
   n<-length(x); xb<-mean(x)
   if (sigma>=0){
      if (side<0){
         tmp<-sigma/sqrt(n)*qnorm(1-alpha)
         a <- -Inf; b <- xb+tmp
      }
      else if (side>0){
         tmp<-sigma/sqrt(n)*qnorm(1-alpha)
         a <- xb-tmp; b <- Inf
      }
      else{
         tmp <- sigma/sqrt(n)*qnorm(1-alpha/2)
         a <- xb-tmp; b <- xb+tmp
      }
      df<-n
   }
   else{
      if (side<0){
         tmp <- sd(x)/sqrt(n)*qt(1-alpha,n-1)
         a <- -Inf; b <- xb+tmp
      }
      else if (side>0){
         tmp <- sd(x)/sqrt(n)*qt(1-alpha,n-1)
         a <- xb-tmp; b <- Inf
      }
      else{
         tmp <- sd(x)/sqrt(n)*qt(1-alpha/2,n-1)
         a <- xb-tmp; b <- xb+tmp
      }
      df<-n-1
   }
   data.frame(mean=xb, df=df, a=a, b=b)
}
