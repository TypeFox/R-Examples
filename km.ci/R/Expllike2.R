"Expllike2" <-
function(lambda,opar)
{
 d<-sum(opar[[2]])
 return(d*log(lambda)-lambda*sum(opar[[1]]))
}

