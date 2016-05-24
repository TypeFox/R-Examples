X.prod<-function(a,b)
{
####   return the cartesian Cross product between two vectors
  z=rep(0, length(a))

  z[1] <- (a[2] * b[3]) - (b[2] * a[3]) 
  z[2] <- (a[3] * b[1]) - (b[3] * a[1]) 
  z[3] <- (a[1] * b[2]) - (b[1] * a[2])
  return(z)

}
