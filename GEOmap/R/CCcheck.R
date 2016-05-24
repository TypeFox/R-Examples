CCcheck<-function(Z)
  {
   #######  make a closed polygon
    x = c(Z$x, Z$x[1])
    y = c(Z$y, Z$y[1])

    n = length(x)
    if(n<=2) return(NULL)


### plot(c(160, 200),c(-85, 85), type='n')
### points(Y)
### lines(Y)
  ###  arrows(Z$x, Z$y,Z$x+2*dx, Z$y+2*dy)

    i = 2:(n-1)
    
    area = sum( ( (x[i]-x[1])*(y[i+1]-y[1]) -    (x[i+1]-x[1])*(y[i]-y[1]) )) 
    
 
    ###  if sign(area) is negative  (<0)direction is clockwise
    ####  if sign(area) is positive (>0) direction is counter-clockwise)
    
    return(sign(area) )

  }
