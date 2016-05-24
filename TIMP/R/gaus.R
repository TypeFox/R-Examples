"gaus" <- 
function(numax, deltanu, nu, nupower=1)
{
  arg <- (2 * (nu - numax))/deltanu
  if(nupower!=1)
    res<-nu^nupower * exp( - log(2) * ((2 * (nu - numax))/deltanu)^2)
  else 
    res<-exp( - log(2) * ((2 * (nu - numax))/deltanu)^2)
  
  res
}

