plot.hnp <-
function(x, cex, pch, colour, lty, lwd, type, xlab, ylab, main, legpos, legcex, ...) {
    if(class(x)!="hnp") stop("Object must be of class 'hnp'")
    if(missing(xlab)) xlab <- "Theoretical quantiles"
    if(missing(ylab)) ylab <- "Residuals"
    if(missing(main)) main <- ""
    if(missing(cex)) cex <- .5
    if(missing(pch)) pch <- 1
    if(missing(colour)) colour <- c(1,1,1,1)
    if(missing(lty)) lty <- c(1,2,1,1)
    if(missing(lwd)) lwd <- c(1,1,1,1)
    if(missing(type)) type <- c("l","l","l","p")
    matplot(x$x, cbind(x$lower, x$median, x$upper, x$residuals),
            col=colour, pch=pch, cex=cex, type=type, lty=lty, lwd=lwd, xlab=xlab, ylab=ylab, main=main, ...)
    if(x$how.many.out) {
      p.out <- round(x$out/x$total*100, 2)
      if(x$print.on) {
        if(missing(legpos)) legpos <- "topleft"
        if(missing(legcex)) legcex <- 1
        legend(legpos, c(paste("Total points:", x$total), paste("Points out of envelope:", x$out, "(", p.out, "%)")), bty="n", cex=legcex)
      } else {
        cat("Total points:", x$total, '\n')
        cat("Points out of envelope:", x$out, "(", p.out, "%)", '\n')
      }
      if(x$paint.out) {
        points(x$out.index[,2], x$out.index[,1], pch=pch, cex=cex, col=x$col.paint.out)
      }
    }
  }
