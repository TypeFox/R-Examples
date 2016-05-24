var_test1<-function(x, sigma2=1, mu=Inf, side=0){
   ## source("p_value.R")
   n<-length(x)
   if (mu<Inf){
      S2<-sum((x-mu)^2)/n; df=n
   }
   else{
      S2<-var(x); df=n-1
   }
   chi2<-df*S2/sigma2;
   P<-p_value(pchisq, chi2, paramet=df, side=side)
   data.frame(var=S2, df=df, chisq2=chi2, P_value=P)
}
