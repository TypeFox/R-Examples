interval_var4<-function(x,y, 
   mu=c(Inf, Inf), side=0, alpha=0.05){ 
   n1<-length(x); n2<-length(y) 
   if (all(mu<Inf)) {
      Sx2<-1/n1*sum((x-mu[1])^2); df1<-n1
      Sy2<-1/n2*sum((y-mu[2])^2); df2<-n2
   }
   else{
      Sx2<-var(x); Sy2<-var(y); df1<-n1-1; df2<-n2-1
   }
   r<-Sx2/Sy2
   if (side<0) {
      a <- 0
      b <- r/qf(alpha,df1,df2)
   }
   else if (side>0) {
      a <- r/qf(1-alpha,df1,df2)
      b <- Inf
   }
   else{
      a<-r/qf(1-alpha/2,df1,df2)
      b<-r/qf(alpha/2,df1,df2)
   }
   data.frame(rate=r, df1=df1, df2=df2, a=a, b=b)
}
