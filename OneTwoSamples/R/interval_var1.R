interval_var1<-function(x,mu=Inf,alpha=0.05){ 
   n<-length(x) 
   if (mu<Inf){
      S2 <- sum((x-mu)^2)/n; df <- n
   }
   else{
      S2 <- var(x); df <- n-1
   }
   a<-df*S2/qchisq(1-alpha/2,df)
   b<-df*S2/qchisq(alpha/2,df)
   data.frame(var=S2, df=df, a=a, b=b)
}
