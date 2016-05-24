"Expllike" <-
function(lambda,times,cens)
{
 d<-sum(cens)
 return(d*log(lambda)-lambda*sum(times))
}

