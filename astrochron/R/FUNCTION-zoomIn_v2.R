### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### zoomIn: explore cross-plot, allow user to zoom into specfied region.
###                                      (SRM: April 16, 2014; June 30, 2014)
#
###########################################################################

zoomIn <- function (dat1,dat2=NULL,ptsize=1,xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,plotype=1,verbose=T)
{
    
    if(verbose) cat("\n----- DYNAMICALLY EXPLORE CROSS-PLOT -----\n")
    dat1 <- data.frame(dat1)
    if(!is.null(dat2)) dat2 <- data.frame(dat2)
    
    if(length(dat1) != 2 && is.null(dat2)) stop("**** TERMINATING: you must input 2 variables")
    if(length(dat1) != 2 && !is.null(dat2)) 
      {
        if(length(dat1) != 1 || length(dat2) != 1) stop("**** TERMINATING: you must input 2 variables")    
      }
    if(length(dat1) == 2 && !is.null(dat2)) stop("**** TERMINATING: you must input 2 variables") 
    
    if(is.null(dat2))
      {
        xx<-dat1[,1]
        yy<-dat1[,2]
      } 
    
    if(!is.null(dat2))
      { 
        xx <- dat1[,1]
        yy <- dat2[,1]
      }  
    
    if(length(xx) < 2) stop("**** TERMINATING: input must have more than one point")
    if(length(xx) != length(yy)) stop("**** TERMINATING: the two variables must have the same number of points")
    
    if(verbose) cat(" * Select region to zoom by clicking\n")
    if(verbose) cat("   Stop by pressing ESC-key (Mac) or STOP button (Windows)\n")
    
    if(is.null(xmin)) xmin=min(xx)
    if(is.null(xmax)) xmax=max(xx)
    if(is.null(ymin)) ymin=min(yy)
    if(is.null(ymax)) ymax=max(yy)
    par(mfrow=c(1,1))
    if (plotype == 1) { plot(xx,yy, main="Click on two spots to zoom (ESC-key or STOP to exit)",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1,cex=ptsize); lines(xx,yy,col="red") }
    if (plotype == 2) { plot(xx,yy, main="Click on two spots to zoom (ESC-key or STOP to exit)",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",cex.axis=1.1,cex.lab=1.1,cex=ptsize) }
    if (plotype == 3) { plot(xx,yy, type="l", main="Click on two spots to zoom (ESC-key or STOP to exit)",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1) }

    doit <- function(x, y=NULL, n=length(x), pch=19, plotype=plotype, cex, ...)
{
    xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
    sel <- rep(FALSE, n); res <- integer(0)
    while(sum(sel) < n) {
        ans <- locator(n=2, type="p",col="blue")
        if(!length(ans)) break
        xmin=min(ans$x)
        xmax=max(ans$x)
        ymin=min(ans$y)
        ymax=max(ans$y)
        if (plotype == 1) { plot(x,y, main="Click on two spots to zoom (ESC-key or STOP to exit)",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1,cex=ptsize); lines(x,y,col="red") }
        if (plotype == 2) { plot(x,y, main="Click on two spots to zoom (ESC-key or STOP to exit)",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",cex.axis=1.1,cex.lab=1.1,cex=ptsize) }
        if (plotype == 3) { plot(x,y, type="l", main="Click on two spots to zoom (ESC-key or STOP to exit)",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1) }
    }
}
   
     doit(xx,yy, cex=ptsize, plotype=plotype)

### END function zoomIn
}