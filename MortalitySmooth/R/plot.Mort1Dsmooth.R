plot.Mort1Dsmooth <-
function(x, 
                              Type=c("logrates", "deaths"), ...){
  object <- x
  Type <- match.arg(Type)
  x <- object$x
  y <- object$y
  Plot <- switch(Type,
                 deaths = 1,
                 logrates = 2)
  if(Plot==1){
    y.hat <- object$fitted.values
    y[object$w==0] <- NA
    matplot(x, cbind(y, y.hat),
            ylab="deaths",
            type=c("p", "l"), col=c(1,2),
            pch=c(1,-1),
            lty=c(0,1),lwd=c(1,2))
  }
  if(Plot==2){
    eta <- log(y) - object$offset
    eta[object$w==0] <- NA
    eta.hat <- object$B %*% object$coef
    matplot(x, cbind(eta, eta.hat),
            ylab="logrates",
            type=c("p", "l"), col=1:2,
            pch=c(1,-1),
            lty=c(0,1),lwd=c(1,2))
    }
}
