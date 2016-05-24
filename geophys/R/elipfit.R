elipfit<-function(ex,ey, PLOT=FALSE, add=TRUE, ...)
  {
    ##   find the least squares ellipse fit to a a set of points
    if(missing(PLOT)) PLOT=FALSE
    if(missing(add)) add=TRUE
    
    mx = mean(ex)
    my = mean(ey)
    kex = ex-mx
    key = ey-my
    
    b = rep(1, length(kex))
    A = cbind(kex^2, kex*key, key^2)
    AtA = t(A) %*% A


    sol = solve(AtA) %*% t(A) %*% b

    LIP = lipper(sol[1], sol[2]/2,sol[3]) 

    if(PLOT)
      {
        theta  = seq(from=0, to=2*pi, length=360)

        phi = LIP[3]
        px = LIP[1]*cos(theta)*cos(phi)-  LIP[2]*sin(theta)*sin(phi)
        py = LIP[1]*cos(theta)*sin(phi) + LIP[2]*sin(theta)*cos(phi)

       if(add==FALSE) plot(ex, ey, type='p', asp=1, ...)
        lines(px+mx, py+my, ...)
      }

    return( LIP )

  }

