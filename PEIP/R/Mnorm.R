Mnorm <-
function(X, k=2)
  {
### the norm function in matlab
##  returns the largest singular value of the matrix or vector
##  here we create a function to duplicate this bizare behavior
    if(missing(k))  k=2

    if(k==1)
      {
        if(is.matrix(X))
          {
            s1 = apply(abs(X), 2, sum)
            return(max(s1))
          }
        else
          {

            s = sum(abs(X))
             return(s)
          }
            
      }

    if(is.matrix(X))
      {
        s = svd(X, 0, 0)
        return(max(s$d))
      }
    else
      {
       s = sqrt(sum(X^2))
       return(s)
      }
    
   

  }
