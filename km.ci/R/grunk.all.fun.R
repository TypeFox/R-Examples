"grunk.all.fun" <-
function(survi,alpha=0.05)
{
   survi <- survi
   
   n.event<- survi$n.event
   n.lost <- survi$n.risk-c(survi$n.risk[-1],0)
   n.cens <- n.lost - n.event
   time <- rep(survi$time,n.lost)

   len <- length(survi$time)
   upper <- numeric(len)
   lower <- numeric(len)
   status <- numeric()
   for(i in 1:length(n.cens))
   {
status <- c(status, rep(0,n.cens[i]),rep(1,n.event[i]))
   }

   for (i in 1:len)
   {
        if(survi$n.event[i]!=survi$n.risk[i])
{
grunk.estimate <- lrt.confints(time,status,survi$time[i],alpha)
        upper[i] <-grunk.estimate$upper
       lower[i] <-grunk.estimate$lower
}
if(survi$n.event[i]==survi$n.risk[i])
{
        upper[i] <-NA
       lower[i] <-NA
}
   }
   survi$upper <- upper
   survi$lower <- lower
   return(survi)

}

