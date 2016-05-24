`CROSSL` <-
function(A1, A2)
{

  if(!is.numeric(A1$x)) A1 = tocartL(A1)
  if(!is.numeric(A2$x)) A2 = tocartL(A2)
    
  x = A1$y*A2$z-A1$z*A2$y
  y = A1$z*A2$x-A1$x*A2$z
  z = A1$x*A2$y-A2$x*A1$y
  a = TOSPHERE(x, y, z)
  return(list(x=x, y=y, z=z, az=a$az, dip=a$dip))
}

