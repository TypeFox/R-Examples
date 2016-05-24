"calcC" <- 
  function (theta, t) 
{
  tfun <- function(t,r) exp(-r*t)
  mapply(tfun, r=theta, MoreArgs=list(t=t)) 
  
}
