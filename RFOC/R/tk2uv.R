tk2uv<-function(T,k)
{
   tau=T*(1-abs(k))
   if(tau >= 0 & k <= 0) {u=tau;v=k}
   if(tau <= 0 & k >= 0) {u=tau;v=k}
   if(tau > 0 & k > 0)
   {
      if(tau <= 4*k)
      {
         u=tau/(1-tau/2);v=k/(1-tau/2)
      }
      if(tau > 4*k)
      {
         u=tau/(1-2*k);v=k/(1-2*k)
      }
   }
   if(tau < 0 & k < 0)
   {
      if(tau >= 4*k)
      {
         u=tau/(1+tau/2);v=k/(1+tau/2)
      }
      if(tau < 4*k)
      {
         u=tau/(1+2*k);v=k/(1+2*k)
      }
   }
   return(list(u=u,v=v))
}
