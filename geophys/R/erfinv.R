`erfinv` <-
  function(x)
  {
    a = (x+1)/2
   w =  qnorm(a)/sqrt(2)
    return( w )
  }

