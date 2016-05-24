var_test2<-function(x, y, mu=c(Inf, Inf), side=0){
   ## source("p_value.R")
   n1<-length(x); n2<-length(y)
   if (all(mu<Inf)){
      Sx2<-sum((x-mu[1])^2)/n1; Sy2<-sum((y-mu[2])^2)/n2
      df1=n1; df2=n2
   }
   else{
      Sx2<-var(x); Sy2<-var(y); df1=n1-1; df2=n2-1
   }
   r<-Sx2/Sy2
   P<-p_value(pf, r, paramet=c(df1, df2), side=side)
   data.frame(rate=r, df1=df1, df2=df2, F=r, P_value=P)
}

