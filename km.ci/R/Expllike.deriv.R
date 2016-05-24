"Expllike.deriv" <-
function(lambda,times,cens)
{
 d<-sum(cens)
 return(d/lambda-sum(times))
}

