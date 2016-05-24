projCircle <- function(x,y,xprev,yprev,circle="row"){
  # projects points of x onto circle if circle = "row"
  ndim <- ncol(x)
  if (circle == "row"){
    meanx <- colMeans(x)
    #x <- x-outer(rep(1,nrow(x)),meanx)
    #y <- y-outer(rep(1,nrow(y)),meanx)
    xunc <- x
    x <- x*outer(rowSums(x^2)^-.5,rep(1,ndim),"*")
    radius <- nrow(x)^-1 * sum(x*xunc)
    x <- radius * x
  } else {
    meany <- colMeans(y)
    #x <- x-outer(rep(1,nrow(x)),meany)
    #y <- y-outer(rep(1,nrow(y)),meany)
    yunc <- y
    y <- y*outer(rowSums(y^2)^-.5,rep(1,ndim),"*") 
    ssOld <- sum((yprev-yunc)^2)
    ssNew <- sum((y-yunc)^2)
    radius <- nrow(y)^-1 * sum(y*yunc)
    y <- radius * y      
  }
  return(list(x=x, y=y))
}
