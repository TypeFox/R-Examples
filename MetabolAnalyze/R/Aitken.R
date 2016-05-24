Aitken <-
function(ll, lla, v, q, epsilon)
{
   la<-0
   tol <- epsilon+1
   if(v > 2)
   {
      av <- (ll[q,v]-ll[q,v-1])/(ll[q,v-1]-ll[q,v-2])  
      lla[q,v]<-ll[q,v-1] + (ll[q,v]-ll[q,v-1])/(1-av)
      la<-lla[q,v]      ## Asymptotic maximum of the data loglikelihood for iteration v.
      tol <- abs(lla[q,v]-lla[q,v-1])
   } #if loop
   list(tol=tol, la=la)
}

