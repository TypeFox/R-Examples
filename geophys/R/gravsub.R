centroid<-function(p)
  {
    n <- length(p$x)
    x1 <- p$x
    i2 <- c(n, 1:(n - 1))
    x2 <- p$x[i2]
    y1 <- p$y
    y2 <- p$y[i2]
    a <- x1 * y2 - x2 * y1
    s <- sum(a) * 3
    if (s == 0) 
      c(mean(x1), mean(y1))
    else c(sum((x1 + x2) * a)/s, sum((y1 + y2) * a)/s)
  }

flipZEE<-function(pol)
  {
    pol$y = -pol$y
    return(pol)
  }

dircheck<-function(pol)
  {
    ####  go around a polygon and
    #### find cross product of successive vectors
    ###  if they are clockwise (Right Handed) they will be positive
    ###  code is vectorized (no loops)
    N = length(pol$x)

    if(N<3) { return(NA) }

    I1 = 1:(N)
    I2 = c(2:(N), 1)
    I3 = c(3:(N), 1, 2)
    
    b = rep(NA, N)
    v1 = cbind(pol$x[I1], pol$y[I1])
    v2 = cbind(pol$x[I2], pol$y[I2])
    v3 = cbind(pol$x[I3], pol$y[I3])
 
    V1x = v3[,1]-v2[,1]
    V2x = v2[,1]-v1[,1]
    V1y = v3[,2]-v2[,2]
     V2y = v2[,2]-v1[,2]
    zp = V1x*V2y - V2x*V1y
    b = sign(zp)
    return(b)
  }


rev2RH<-function(pol)
  {
    pol = list(x=rev(pol$x), y= rev(pol$y))
    return(pol)
    
  }
