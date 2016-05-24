
OverlapCurve <- function(object, plot = TRUE, xlim , ylim , xlab , ylab , main , ...)  {
   # object - an object of class "rvals"
   
   if(missing(xlim)) {
       xlim = range(object$aux$alpha.grid)
   }
   if(missing(ylim)) {
       ylim = c(0,1)
   }
   if(missing(xlab)) {
       xlab = "proportion selected"
   }
   if(missing(ylab)) {
       ylab = "overlap"
   }
   if(missing(main)) {
       main = "Estimated overlap curve"
   }
   
   ovc <- object$aux$alpha.grid
   ngrid <- length(ovc)
   for(j in 1:ngrid) {
        thresh <- (object$aux$unsorted$RValue <= object$aux$alpha.grid[j])
        object$aux$alpha.grid
        ovc[j] <- mean(object$aux$V[,j]*thresh)
   }
   ans <- approxfun(object$aux$alpha.grid,ovc)
   if(plot)  {
      plot(object$aux$alpha.grid,ovc,type="n",xlim=xlim,ylim=ylim,xlab=xlab,
           ylab=ylab,main=main, ...)
      lines(object$aux$alpha.grid,ovc,lwd=2)
      lines(object$aux$alpha.grid,object$aux$alpha.grid,lwd=2,lty=2)
   
      invisible(ans)
   }
   else {
      return(ans)
   }
}
