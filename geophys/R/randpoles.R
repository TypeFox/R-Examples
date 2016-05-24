randpoles<-function(az ,   iang, alphadeg, opt="unif", BALL.radius = 1,N = 10,  add=TRUE, ...)
  {
    if(missing(BALL.radius)) BALL.radius = 1
    if(missing(add))  add=TRUE
    if(missing(N))  N = 10
    
    DEG2RAD = pi/180
##########   net(); addsmallcirc(30, 40, 10)
    if(missing(opt)) { opt="unif" }

    phi =  runif(N, 0, 2*pi)
    theta = runif(N, 0, pi*alphadeg/180)
    
    
    if(opt=="unif")
      {
        phi =  runif(N, 0, 2*pi)
        theta = runif(N, 0, pi*alphadeg/180)
      }

    
    if(opt=="norm")
      {
        phi =  rnorm(N, mean=0, sd=0.75*(2*pi))
        theta = rnorm(N, mean=0, sd=0.75*(pi*alphadeg/180))
      }

    x = BALL.radius*sin(theta)*cos(phi)
    y =  BALL.radius*sin(theta)*sin(phi)
    z = BALL.radius*cos(theta)
    
    D = cbind(x,y,z)
    

    ry =  RFOC::roty3(iang)
    
    rz =  RFOC::rotz3(az)

    Rmat = ry %*% rz
    g = D %*% Rmat
    
    r2 = sqrt(g[,1]^2+g[,2]^2+g[,3]^2)

    phi2 = atan2(g[,2], g[,1])
    
    theta2 = acos(g[,3]/r2)


    raz =  phi2*180/pi
    rdip = theta2*180/pi

    
    Sc = qpoint(phi2*180/pi,  theta2*180/pi, PLOT=FALSE)
    
    diss = sqrt((Sc$x[1:(N-1)]-Sc$x[2:(N)])^2+(Sc$y[1:(N-1)]-Sc$y[2:(N)])^2)
    
    ww = which(diss>0.9*BALL.radius)
    if(add==TRUE)
      {
        
        points(Sc$x, Sc$y, ...)
        
      }
    invisible(list(az=raz, dip=rdip, x=Sc$x, y=Sc$y ))

  }
