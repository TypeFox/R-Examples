p_value<-function(cdf, x, paramet=numeric(0), side=0){
   n<-length(paramet)
   P<-switch(n+1,
      cdf(x), 
      cdf(x, paramet), 
      cdf(x, paramet[1], paramet[2]),
      cdf(x, paramet[1], paramet[2], paramet[3])
   )
   if (side<0)         P
   else if (side>0)  1-P
   else 
      if (P<1/2)     2*P 
      else           2*(1-P)
}