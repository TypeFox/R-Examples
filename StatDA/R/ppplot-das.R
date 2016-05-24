ppplot.das <- function(x, pdist=pnorm, 
                      xlab=NULL, ylab="Probability", line=TRUE, lwd=2, pch=3, 
                      cex=0.7, cex.lab=1, ...)
{
    DOTARGS <- as.list(substitute(list(...)))[-1]
    DOTARGS <- paste(names(DOTARGS), DOTARGS, sep="=", collapse=", ")

#    xlab=deparse(substitute(x))

    PNAME <- deparse(substitute(pdist))
    pdist <- match.fun(pdist)
    
    y <- ppoints(length(x))

    pprobs <- pdist(sort(x), mean(x), sd(x), ...)
      
    plot(y, pprobs, axes=FALSE, type="n", xlab=xlab, ylab=ylab,
         xlim=c(0,1), ylim=c(0,1), cex.lab=cex.lab)
    box()
   
    probs <- seq(0,1,by=0.1) 
    axis(1)
    axis(2)

    points(y, pprobs, pch=pch, cex=cex)

#    PTEXT <- paste("Probabilities: ", PNAME, sep="")
#    if(nchar(DOTARGS))
#        PTEXT <- paste(PTEXT, DOTARGS, sep=", ")
#    mtext(PTEXT, side=1, line=3, adj=1)
    
    if(line){
        yl <- quantile(y, c(0.25, 0.75))
        pl <- quantile(pprobs, c(0.25, 0.75))
        slope <- diff(pl)/diff(yl)
        int <- pl[1] - slope * yl[1]
        abline(int, slope, col=1, lwd=lwd)
    }
}
