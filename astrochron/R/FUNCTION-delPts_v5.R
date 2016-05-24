### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### delPts: interactively delete points in plot (SRM: June 10-11, 2013
###                                           June 13, 2013; June 19, 2013;
###                                           July 27, 2013; April 3, 2014;
###                                           June 15,  2015)
####
###########################################################################

delPts <- function (dat,del=NULL,ptsize=1,xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,plotype=1)
{
   
    cat("\n----- INTERACTIVELY IDENTIFY AND DELETE POINTS IN PLOT -----\n")
    dat <- data.frame(dat)
    if(length(dat) != 2) {stop("**** TERMINATING: input must have two columns")}
    if(nrow(dat) < 2) {stop("**** TERMINATING: input must have more than one point")}
    
    ipts <- length(dat[,1]) 
    cat(" * Number of data points=", ipts,"\n")
    
# skip interactive plot and delete points if del=!NULL
    if(!is.null(del))  pts=del

# otherwise proceed with interactive plot
    if(is.null(del))
     { 
       cat(" * Select points by clicking\n")
       cat("   Stop by pressing ESC-key (Mac) or STOP button (Windows)\n")
    
       if(is.null(xmin)) xmin=min(dat[,1])
       if(is.null(xmax)) xmax=max(dat[,1])
       if(is.null(ymin)) ymin=min(dat[,2])
       if(is.null(ymax)) ymax=max(dat[,2])
       par(mfrow=c(1,1))
       if (plotype == 1) { plot(dat, main="Select data points for deletion (press ESC-key or STOP to exit)",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1,cex=ptsize); lines(dat,col="red") }
       if (plotype == 2) { plot(dat, main="Select data points for deletion (press ESC-key or STOP to exit)",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",cex.axis=1.1,cex.lab=1.1,cex=ptsize) }
       if (plotype == 3) { plot(dat, type="l", main="Select data points for deletion (press ESC-key or STOP to exit)",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1) }

## this script modified from '?identify' in R
identifyPch <- function(x, y=NULL, n=length(x), pch=19, cex, ...)
{
    xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
    sel <- rep(FALSE, length(x)); res <- integer(0)
    while(sum(sel) < n) {
# note: plot must be set to FALSE, as numbers plotted are from culled data set!
        ans <- identify(x[!sel], y[!sel], n=1, plot=F, ...)
#        ans <- identify(x[!sel], y[!sel], n=n, plot=T, ...)
        if(!length(ans)) break
        ans <- which(!sel)[ans]
        points(x[ans], y[ans], pch = pch, cex = cex, col="blue")
        sel[ans] <- TRUE
        res <- c(res, ans)
    }
    res
}

       pts <- identifyPch(dat[,1],dat[,2], cex=ptsize)
   }
 
    cat("\nSELECTED DATA POINTS FOR DELETION:\n")
    print(dat[pts,])
    cat("\n")

    out <- dat
    out[pts,] <- NA
    out <- data.frame(subset(out, !(out[, 2] == "NA")))

    newpts=length(out[,1])
    cat(" * Number of data points following deletion=",newpts,"\n")

    par(mfrow=c(2,1))
    if (plotype == 1) { plot(dat, main="Original Data Series (deleted points in blue)",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1,cex=ptsize/1.5,col="gray"); lines(dat,col="red") }
    if (plotype == 2) { plot(dat, main="Original Data Series (deleted points in blue)",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",cex.axis=1.1,cex.lab=1.1,cex=ptsize/1.5,col="gray") }
    if (plotype == 3) { plot(dat, type="l", main="Original Data Series (deleted points in blue)",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1,col="gray") }
    points(dat[pts,1],dat[pts,2],col="blue",pch=19,cex=ptsize/1.5)

    if (plotype == 1) { plot(out, main="Edited Data Series",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1,cex=ptsize/1.5,col="gray"); lines(out,col="red") }
    if (plotype == 2) { plot(out, main="Edited Data Series",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",cex.axis=1.1,cex.lab=1.1,cex=ptsize/1.5,col="gray") }
    if (plotype == 3) { plot(out, type="l", main="Edited Data Series",xlim=c(xmin,xmax),ylim=c(ymin,ymax),bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1,col="gray") }

    return(out)
### END function delPts
}
