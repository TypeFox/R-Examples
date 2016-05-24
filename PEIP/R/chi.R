chi <-
function(x,n)
  {
    f = ( 1 / (2^(n/2)*gamma(n/2)) )* (x^( (n/2) -1))* exp(-x/2)
    return(f)
  }
