epicurve.table <- function(x,
                           width = 1,
                           space = 0,
                           tick = TRUE,
                           tick.offset = 0.5,
                           segments = FALSE,
                           ...){
  xvals <- barplot(x, width=width, space=space, ...)
  if(tick) {axis(1, at=c(0, xvals + tick.offset), labels=FALSE, tick=TRUE)}
  if(segments){
    xx <- xvals-(width/2)
    y2 <- apply(x,2,sum)
    xy2 <- cbind(xx,y2)
    y0 <- cbind(xy2[1,1],0:xy2[1,2])
    z0 <- cbind(y0, y0[,1]+width, y0[,2])
    for(i in 2:nrow(xy2)){
      yy <- cbind(xy2[i,1],0:xy2[i,2])
      z <- cbind(yy, yy[,1]+width, yy[,2])
      z2 <- rbind(z0,z)
      z0 <- z2
    }
    segments(z0[,1],z0[,2],z0[,3],z0[,4])
  }
  invisible(xvals)
}

