`xprod` <-
function(A1,A2)
  {
############  cross product with 2 vectors input
    x = A1[2]*A2[3]-A1[3]*A2[2]
    y = A1[3]*A2[1]-A1[1]*A2[3]
    z = A1[1]*A2[2]-A2[1]*A1[2]
   ###  return a vector
    return(c(x,y,z))
  }

