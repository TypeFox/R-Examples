`gpoly` <-
function(x)
  {
    
    e = x;  
    n = length(e);
    c = c(1, rep(0,n));
    for(j in 1:n)
      {
        c[2:(j+1)] = c[2:(j+1)] - e[j]*c[1:j];
      }
    return(Re(c))
  }

