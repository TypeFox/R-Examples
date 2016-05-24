philpotts.erf <-
function(x)
  {
    ## copied from phillpots
    
    t = 1/(1+0.47047*x)
    a1 = 0.34802
    a2 = -0.09587
    a3 = 0.74785

    y = 1 - (a1*t + a2*t^2 + a3*t^3)* exp(-(x^2))

    return(y)

  }

