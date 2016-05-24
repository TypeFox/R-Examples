`AXpoint` <-
function(UP=TRUE, col=2, n=1)
{
                   #   to plot points on an equal area stero net:
  if(missing(UP)) { UP = TRUE}
  if(missing(col)) {col='red' }
  if(missing(n)) { n = NULL }
  
  deg2rad = pi/180
  rad2deg = 180/pi
  
  if(is.null(n))
     { p = locator(type='p', col=col) }
  else
    { p = locator(n=n, type='p', col=col)  }


  r = sqrt(p$x^2+p$y^2)

  goodones = r<=1
  
  iang = rad2deg*2*asin(r[goodones]/sqrt(2))
  phiang = rad2deg*( pi/2 - atan2(p$y[goodones],p$x[goodones]))
  
  if(UP==TRUE)
    {
      iang = 180-iang
    }
  a = RSEIS::TOCART(phiang, iang)
  return(list(az=phiang, dip=iang, x=a$x, y=a$y, z=a$z, gx=p$x, gy=p$y))
}

