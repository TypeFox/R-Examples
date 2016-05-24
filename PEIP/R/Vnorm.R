Vnorm <-
function(X)
  {
### the norm function in matlab
##  returns the largest singular value of the matrix or vector
##  here we create a function to duplicate this bizare behavior    
    
    return( sqrt(sum(X^2)) )

  }
