trackSpeed <- function(coords, time) {
    dt <- diff(as.numeric(time))
    d2 <- diff(as.matrix(coords))^2
    speed <- d2[,1]
    ndim <- ncol(d2)
    if(ndim > 1)
      speed <- d2 %*% rep(1,ndim)
    else
      speed <- d2
    sqrt(as.vector(speed))/dt
}
    
