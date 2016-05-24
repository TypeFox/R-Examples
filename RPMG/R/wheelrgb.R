`wheelrgb` <-
function(wloc, v, RY)
  {
    ###########  function used in colwheel
    
    upp = wloc$y>RY[2]
    
    ay = wloc$y[upp] - diff(RY)
    ax = wloc$x[upp]
    
    mcol = list(x=ax, y = ay)
    
    mARE = sqrt((mcol$x)^2+(mcol$y)^2)
    
    mARE[mARE>1] = 1
    
    mANG = atan2(mcol$y,mcol$x)
    
    mpANG = (mANG+pi)/(2*pi)
    
    colsupp = hsv( mpANG, v, mARE)
    
    
    ay = wloc$y[!upp] 
    ax = wloc$x[!upp]
    
    mcol = list(x=ax, y = ay)
    
    mARE = sqrt((mcol$x)^2+(mcol$y)^2)
    
    mARE[mARE>1] = 1
    
    mANG = atan2(mcol$y,mcol$x)
    
    mpANG = (mANG+pi)/(2*pi)
    
    colsdwn = hsv( mpANG, mARE, v)
    
    
    zcols = rep(NA, length(wloc$x))
    zcols[upp] = colsupp
    zcols[!upp] = colsdwn
    
    return(zcols)
    
  }

