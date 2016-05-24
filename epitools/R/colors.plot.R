"colors.plot" <-
  function(locator = FALSE, cex.axis = 0.7){
    xx <- rep(1:30,22)
    yy <- rep(1:22,rep(30,22))
    yyy <- matrix(yy,ncol=22)
    cm <- colors.matrix()
    matplot(xx[1:30], yyy, pch=15, type="n", axes=FALSE,
            xlab="colors.matrix[row, ]",
            ylab="colors.matrix[ , col]",
            main ="Matrix plot of 'colors()' function in R.
            Use coordinates to identify color name.")
    title(sub = "Source: www.epitools.net", cex.sub = 0.7)
    points(xx,yy, type="p", pch=15, cex=2, col = c(colors(),NA,NA,NA))
    axis(1, at=c(0:30 + 0.5), labels=FALSE, tick=TRUE)
    axis(1, at=1:30, labels=1:30, cex.axis=cex.axis, tick=FALSE)
    axis(2, at=c(0:22 + 0.5), labels=FALSE, tick=TRUE)
    axis(2, at=1:22, labels=1:22, cex.axis=cex.axis, tick=FALSE, las=1)
    if(locator==TRUE){
      lxy <- locator()
      xy <- round(data.frame(lxy))
      xym <- as.matrix(xy)
      located <- data.frame(xy, color.names = cm[xym])
      return(located)
    } else invisible(cm)
  }
