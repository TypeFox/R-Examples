plot.abcrf <- function(x, obs=NULL, n.var=20, main="", pdf=FALSE, ...)
{
	old.par <- par(no.readonly = TRUE)
	if (length(x$model.rf$importance)<20) n.var <- length(x$model.rf$importance)
	modindex <- x$model.rf$y
	sumsta <- x$sumsta
 	if (x$lda) {
 	  if (pdf) { 
 	    pdf("graph_varImpPlot.pdf")
		  varImpPlot(x$model.rf, n.var=n.var, main=main, ...)
		  dev.off()
 	  }
 	  varImpPlot(x$model.rf, n.var=n.var, main=main, ...)
		nmod <- nlevels(modindex)
		nstat <- ncol(sumsta)-(nmod-1)
		projections <- sumsta[,(nstat+1):(nstat+nmod-1)]
		if  (!is.null(obs)) projobs <- predict(x$model.lda,obs)$x
		coloris <- rainbow(nmod)
		colo <- coloris[modindex]
		readline("Press <ENTER> to Continue")
    if (nmod > 2) {
      if (pdf)
        {
        pdf("graph_lda.pdf")
        plot(projections[,1:2], col=colo, pch=3)
        legend("topleft", legend = as.character(levels(modindex)), col = coloris, 
               pch = 15, bty = "o", pt.cex = 2, cex = .8, horiz = TRUE, 
               inset = c(.01, .01), title = "Models", bg = "white")
	  	  if  (!is.null(obs)) points(projobs[1],projobs[2],pch="*",cex=5.3)
	  	  dev.off()
		    }
      plot(projections[,1:2], col=colo, pch=3)
      legend("topleft", legend = as.character(levels(modindex)), col = coloris, 
             pch = 15, bty = "o", pt.cex = 2, cex = .8, horiz = TRUE, 
             inset = c(.01, .01), title = "Models", bg = "white")
      if  (!is.null(obs)) points(projobs[1],projobs[2],pch="*",cex=5.3)
    } else {
      l1 <- levels(modindex)[1]
      l2 <- levels(modindex)[2]
      d1 <- density(projections[modindex == l1])
      d2 <- density(projections[modindex == l2])
      coloris <- c("blue", "orange")
      xrange <- range(c(d1$x, d2$x))
      yrange <- c(0, 1.2*max(c(d1$y, d2$y)))
      if (pdf)
        {
        pdf("graph_lda.pdf")
        plot(d1, xlim = xrange, ylim = yrange,
             col=coloris[1], main="", xlab="")
        lines(d2, col=coloris[2])
        legend("topleft", legend = as.character(levels(modindex)), col = coloris, 
                cex = .8, horiz = TRUE, lty=1, bty="o",
               inset = c(.01, .01), title = "Models", bg = "white")
      	if  (!is.null(obs)) abline(v=projobs)
        dev.off()
      }
      plot(d1, xlim = xrange, ylim = yrange,
           col=coloris[1], main="", xlab="")
      lines(d2, col=coloris[2])
      legend("topleft", legend = as.character(levels(modindex)), col = coloris, 
              cex = .8, horiz = TRUE, lty=1, bty="o",
             inset = c(.01, .01), title = "Models", bg = "white")
      if  (!is.null(obs)) abline(v=projobs)
    }
	} else {
	  if (pdf)
	    {
	    pdf("graph_varImpPlot.pdf")
	    varImpPlot(x$model.rf, n.var=n.var, ...)
	    dev.off()
	    }
		varImpPlot(x$model.rf, n.var=n.var, main=main, ...)
	}
	par(old.par)
}