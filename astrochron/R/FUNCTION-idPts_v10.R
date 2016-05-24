### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### idPts: identify points in plot (SRM: Oct. 29, 2012; Nov. 23, 2012;
###                                      May 20, 2013; June 10-11, 2013;
###                                      June 13, 2013; July 27, 2013;
###                                      April 10, 2014; April 21, 2014;
###                                      April 23, 2014; Nov. 19, 2014;
###                                      February 4, 2015; June 1, 2015)
###########################################################################


idPts <- function (dat1,dat2=NULL,ptsize=1,xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,logx=F,logy=F,plotype=1,annotate=1,output=1,verbose=T)
{    
    if(verbose) cat("\n----- INTERACTIVELY IDENTIFY POINTS IN PLOT -----\n")
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

    if(verbose) cat(" * Select points by clicking\n")
    if(verbose) cat("   Stop by pressing ESC-key (Mac) or STOP button (Windows)\n")

    if(is.null(xmin)) xmin=min(xx)
    if(is.null(xmax)) xmax=max(xx)
    if(is.null(ymin)) ymin=min(yy)
    if(is.null(ymax)) ymax=max(yy)
   
    par(mfrow=c(1,1))
    if(logx && logy) setlog="xy"
    if(logx && !logy) setlog="x"
    if(!logx && logy) setlog="y"
    if(!logx && !logy) setlog=""
    if (plotype == 1) { plot(xx,yy, main="Select data points (press ESC-key or STOP to exit)",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1,cex=ptsize,xlab="Location",ylab="Value",log=setlog); lines(xx,yy,col="red") }
    if (plotype == 2) { plot(xx,yy, main="Select data points (press ESC-key or STOP to exit)",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",cex.axis=1.1,cex.lab=1.1,cex=ptsize,log=setlog) }
    if (plotype == 3) { plot(xx,yy, type="l", main="Select data points (press ESC-key or STOP to exit)",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1,log=setlog) }

## this script modified from '?identify' in R
identifyPch <- function(x, y=NULL, n=length(x), pch=19, cex, ...)
{
    xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
    sel <- rep(FALSE, length(x)); res <- integer(0)
    while(sum(sel) < n) { 
        ans <- identify(x[!sel], y[!sel], n=1, plot=F, ...)
        if(!length(ans)) break
        ans <- which(!sel)[ans]
        points(x[ans], y[ans], pch = pch, cex = cex, col="blue")
        sel[ans] <- TRUE
        res <- c(res, ans)
        if(annotate==1) text(x[ans],y[ans],round(x[ans],digits=2),col = "blue", pos = 3, offset = 1.2, font=2, xpd = TRUE)
        if(annotate==1) text(x[ans],y[ans],round(y[ans],digits=2),col = "blue", pos = 3, offset = 0.5, font=2, xpd = TRUE)
        if(annotate==2) text(x[ans],y[ans],round(x[ans],digits=2),col = "blue", pos = 1, offset = 1.2, font=2, xpd = TRUE)
        if(annotate==2) text(x[ans],y[ans],round(y[ans],digits=2),col = "blue", pos = 1, offset = 0.5, font=2, xpd = TRUE)
        if(verbose) cat(ans,x[ans],y[ans],"\n")
    }
    res
}
   
    if(verbose) cat("\nSELECTED DATA POINTS:\n")
    pts <- identifyPch(xx,yy, cex=ptsize)

    if(output==1) 
     {
       out <- data.frame(cbind(xx[pts],yy[pts]))
       colnames(out) <- c("x","y")
       return(out)
     }  
        
    if(output==2) 
     {
       out <- data.frame(cbind(pts,xx[pts],yy[pts]))
       colnames(out) <- c("index","x","y")
       return(out)
     }  

### END function idPts
}