int.bspline2 <- function(penden.env,Y) {
  if(is.matrix(Y)) len.Y <- dim(Y)[1]
  if(is.vector(Y)) len.Y <- length(Y)
  x.help <- get("x.help",penden.env)
  cal.int.base <- get("int.base",penden.env)
  int.base <- matrix(0,len.Y,dim(cal.int.base)[2])
  for(j in 1:dim(cal.int.base)[2]) {
    int.base[,j] <- spline(x=x.help,y=cal.int.base[,j],xout=Y)$y
  }
  return(int.base)
}
