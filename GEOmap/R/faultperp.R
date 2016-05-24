`faultperp` <-
function(x,y, N=20,  endtol=.1, h=1, col='black')
  {
    if(missing(N)) { N = 1 }
    if(missing(h)) { h = 1 }
     if(missing(col)) { col='black' }
      if(missing(endtol)) { endtol=.1 }

    g = PointsAlong(x, y, N=N,  endtol=endtol)
    segments(g$x, g$y, g$x+h*(-g$rot$sn)  , g$y+h*(g$rot$cs) , col=col)

    invisible(g)
   
  }

