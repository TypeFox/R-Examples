# panel function for use with coplot (J. Fox)

# last modified 2 April 2009

panel.car <- function(x, y, col, pch, cex=1, span=.5, lwd=2,
    reg.line=lm, lowess.line=TRUE,...){
    points(x, y, col=col, pch=pch, cex=cex)
    if (is.function(reg.line)) regLine(reg.line(y ~ x), 
        lty=2, lwd=lwd, col=col, ...)
    if (lowess.line) lines(lowess(na.omit(as.data.frame(cbind(x, y))), f=span), 
        col=col, lwd=lwd, ...)
    }
