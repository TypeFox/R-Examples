`addsmallcirc` <-
function(az ,   iang, alphadeg, BALL.radius = 1,N = 100,  add=TRUE, ...)
  {
    if(missing(BALL.radius)) BALL.radius = 1
    if(missing(add))  add=TRUE
    if(missing(N))  N = 100
    
      
    DEG2RAD = pi/180
##########   net(); addsmallcirc(30, 40, 40)

    
#######   lat lon of point
#######   radius in degrees
    
   ###   alphadeg = 180*rad /(pi*R.MAPK)


    phi = seq(from=0, to=2*pi, length=N)
    theta = pi*alphadeg/180


    x = BALL.radius*sin(theta)*cos(phi)
    y =  BALL.radius*sin(theta)*sin(phi)
    z = rep(BALL.radius*cos(theta), length(phi))

    D = cbind(x,y,z)


    ry = roty3(iang)

    rz = rotz3(az)

    Rmat = ry %*% rz
    g = D %*% Rmat

    r2 = sqrt(g[,1]^2+g[,2]^2+g[,3]^2)

    phi2 = atan2(g[,2], g[,1])

    theta2 = acos(g[,3]/r2)

Sc = qpoint(phi2*180/pi,  theta2*180/pi, PLOT=FALSE)

diss = sqrt((Sc$x[1:(N-1)]-Sc$x[2:(N)])^2+(Sc$y[1:(N-1)]-Sc$y[2:(N)])^2)

ww = which(diss>0.9*BALL.radius)

    if(length(ww)>0)
      {
        Sc$x[ww] = NA
        Sc$y[ww] = NA

 }
    
    ####Sc = qpoint(nlon, nlat)
####  net()
   #### lines(Sc$x, Sc$y)
    if(add==TRUE)
      {
        
        lines(Sc$x, Sc$y, ...)
        
      }
    invisible(list(x=Sc$x, y=Sc$y ))

  }

