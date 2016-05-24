qpplot.das <- function(x, qdist=qnorm, probs=NULL, logx=FALSE, cex.lab=1,
                 xlab=NULL, ylab="Probability [%]", line=TRUE, lwd=2, pch=3, 
                 logfinetick=c(10),logfinelab=c(10),cex=0.7,xlim=NULL,ylim=NULL,
                 gridy=TRUE, add.plot=FALSE,col=1, ...)
{
#
# logx=TRUE ... use log-scale on x-axis
# logfinetick ... how fine are the tick marks on log-scale on x-axis 
# logfinelab ... how fine are the labels on log-scale on x-axis 
# gridy ... if grid along y-axis should be drawn
# add.plot ... if TRUE the new plot is added to an old one
#
    DOTARGS <- as.list(substitute(list(...)))[-1]
    DOTARGS <- paste(names(DOTARGS), DOTARGS, sep="=", collapse=", ")

#    xlab=deparse(substitute(x))

    x <- sort(x)
    QNAME <- deparse(substitute(qdist))
    qdist <- match.fun(qdist)
    
    y <- qdist(ppoints(length(x)), ...)

    if(is.null(probs)){
        probs <- c(.01, .05, seq(.1,.9, by=.1), .95, .99)
        if(length(x)>=500)
            probs <- c(0.001, probs, .999)
    }

    qprobs <- qdist(probs, ...)

if (!add.plot){
    if (is.null(ylim))
      plot(x, y, axes=FALSE, type="n", ylim=range(c(y,qprobs)),
         xlab=xlab, ylab=ylab, cex.lab=cex.lab, xlim=xlim)
    else
    plot(x, y, axes=FALSE, type="n", ylim=ylim,
         xlab=xlab, ylab=ylab, cex.lab=cex.lab, xlim=xlim)
    box()
}   
    if (gridy){ 
        abline(h=qprobs, lty=3, col=gray(0.5))
    }
    if (logx){
      axis(1,at=log10(alog<-sort(c((10^(-50:50))%*%t(logfinelab)))),labels=alog)
      abline(v=log10(sort(c((10^(-50:50))%*%t(logfinetick)))),lty=3,col=gray(0.5)) 
    }
    else {
      axis(1)
      abline(v=axTicks(1),lty=3,col=gray(0.5)) 
    }
    axis(2, at=qprobs, labels=100*probs)

    points(x, y, pch=pch, cex=cex, col=col)

#    QTEXT <- paste("Quantile: ", QNAME, sep="")
#    if(nchar(DOTARGS))
#        QTEXT <- paste(QTEXT, DOTARGS, sep=", ")
#    mtext(QTEXT, side=1, line=3, adj=1)
    
    if(line){
        xl <- quantile(x, c(0.25, 0.75))
        yl <- qdist(c(0.25, 0.75), ...)
        slope <- diff(yl)/diff(xl)
        int <- yl[1] - slope * xl[1]
        abline(int, slope, col=1, lwd=lwd)
    }
}
