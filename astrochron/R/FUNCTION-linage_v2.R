### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### linage: tune stratigraphic series to astronomical target using interface
###          similar to Analyseries 'Linage' routine (SRM: June 10-13, 2013;
###                                                        June 29-30, 2015)
###
###########################################################################

linage <- function (dat,target,extrapolate=F,xmin=NULL,xmax=NULL,tmin=NULL,tmax=NULL,size=1,plotype=1,output=1,genplot=T)
{
    
    cat("\n----- INTERACTIVELY TUNE STRATIGRAPHIC SERIES TO ASTRONOMICAL TARGET -----\n")
    cat("\n                *****  Select points by clicking   *****\n")
    cat("                *****  Proceed from left to right  *****\n")
    cat("                *****  Stop by pressing ESC-key    *****\n")

    dat <- data.frame(dat)
    if(length(dat) != 2) {stop("**** TERMINATING: stratigraphic series must have two columns")}
    if(nrow(dat) < 2) {stop("**** TERMINATING: stratigraphic series must have more than one point")}
    target <- data.frame(target)
    if(length(target) != 2) {stop("**** TERMINATING: target must have two columns")}
    if(nrow(target) < 2) {stop("**** TERMINATING: target must have more than one point")}

# error checking
    nptsDat <- length(dat[,1]) 
    cat("\n * Number of points in stratigraphic series=", nptsDat,"\n") 
# determine the current direction of the data series   
    dtDat <- dat[2,1]-dat[1,1]
    nptsTarg <- length(target[,1])
    cat(" * Number of points in astronomical target series=", nptsTarg,"\n") 
# determine the current direction of the data series   
    dtTarg <- target[2,1]-target[1,1]

    if(dtDat<0)
     {
       cat("**** Warning: Stratigraphic series is not in increasing order. Sorting now to fix this.\n")
       dat <- dat[order(dat[1],na.last=NA,decreasing=F),]
     }  
    
     if(dtTarg<0)
     {
       cat("**** Warning: Astronomical target series is not in increasing order. Sorting now to fix this.\n")
       target <- target[order(target[1],na.last=NA,decreasing=F),]
     }  
    
    cat("\nSELECTED AGE CORRELATION POINTS (depth/height, time):\n")
    
    if(is.null(xmin)) xmin=min(dat[,1])
    if(is.null(xmax)) xmax=max(dat[,1])
    if(is.null(tmin)) tmin=min(target[,1])
    if(is.null(tmax)) tmax=max(target[,1])


## this script modified from '?identify' in R
identifyPch <- function(x, y=NULL, n=length(x), pch=3, ...)
{
    xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
    sel <- rep(FALSE, length(x)); res <- integer(0)
    while(sum(sel) < n) {
        ans <- identify(x[!sel], y[!sel], n=n, plot=T, ...)
        if(!length(ans)) break
        ans <- which(!sel)[ans]
        points(x[ans], y[ans], pch = pch, col="blue")
        sel[ans] <- TRUE
        res <- c(res, ans)
    }
    res
}

    par(mfrow=c(2,1))

    itick = 0
    ijtick = 0
    test = rep( FALSE,length(dat[,1]) )

# start while loop    
    while( sum(test) < length(dat[,1]) )
     {
# select tuning point from stratigraphic series
       if (plotype == 1) {plot(dat, main="Select stratigraphic correlation points",xlim=c(xmin,xmax),bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1,cex=0.5); lines(dat,col="red")}
       if (plotype == 2) {plot(dat, main="Select stratigraphic correlation points",xlim=c(xmin,xmax),bty="n",cex.axis=1.1,cex.lab=1.1,cex=0.5)}
       if (plotype == 3) {plot(dat, type="l", main="Select stratigraphic correlation points",xlim=c(xmin,xmax),bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1,cex=0.5)}
       if (itick == 0)
         {
           datpts <- identifyPch(dat[,1],dat[,2],n=1)
           ndatpts = 1
         }
       if (itick > 0)
         {
           odatpts=ndatpts
           points(dat[datpts,1],dat[datpts,2],pch=19,col="blue",cex=1*size)
           text(dat[datpts,1],dat[datpts,2],labels=1:itick,cex=0.5*size,col="white",font=2)
           datpts <- append(datpts,identifyPch(dat[,1],dat[,2],n=1))
           ndatpts = length(datpts)    
           if(odatpts == ndatpts) break    
         }
       itick = itick + 1  

# select tuning point from astronomical target  
       if (plotype == 1) {plot(target, main="Select astronomical target correlation points",xlim=c(tmin,tmax),bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1,cex=0.5); lines(target,col="red")}
       if (plotype == 2) {plot(target, main="Select astronomical target correlation points",xlim=c(tmin,tmax),bty="n",cex.axis=1.1,cex.lab=1.1,cex=0.5)}
       if (plotype == 3) {plot(target, type="l", main="Select astronomical target correlation points",xlim=c(tmin,tmax),bty="n",lwd=2,cex.axis=1.1,cex.lab=1.1,cex=0.5)}
       if (ijtick == 0) targetpts <- identifyPch(target[,1],target[,2],n=1)
       if (ijtick > 0)
         {
           points(target[targetpts,1],target[targetpts,2],pch=19,col="blue",cex=1*size)
           text(target[targetpts,1],target[targetpts,2],labels=1:ijtick,cex=0.5*size,col="white",font=2)
           targetpts <- append(targetpts,identifyPch(target[,1],target[,2],n=1))
          } 
       ijtick = ijtick + 1
       
# output correlation points (depth/height vs. time) to screen
       cat(dat[datpts[itick],1],target[targetpts[ijtick],1],"\n")

# end while loop
     }   

# now tune          
    cat("\n * Tuning via piecewise linear interpolation of age control points.\n")
    controlPts = data.frame( cbind(dat[datpts,1],target[targetpts,1]) )
# sort to ensure increasing order
    controlPts <- controlPts[order(controlPts[1],na.last=NA,decreasing=F),]
    colnames(controlPts)[1] = 'Depth/Height'
    colnames(controlPts)[2] = 'Time' 
# note: here we will force tune.R to extrapolate. if extrapolation not desired, we will remove
#   points later. this will allow us to determine the amount necessary to shift control points
#   for plotting purposes
    tuned = tune(dat,controlPts,genplot=F,extrapolate=T,verbose=F)  
    colnames(tuned)[1] = colnames(target[1])
    colnames(tuned)[2] = colnames(dat[2])

    if(extrapolate) remov = 0
    if(!extrapolate) 
     {
# determine number of points removed from start of series, save for plotting purposes
       remov = nrow(subset(tuned, (tuned[1] < target[targetpts[which.min(targetpts)],1]))[1] )
# now remove points from both start and end if needed
       tuned = subset(tuned, (tuned[1] >= target[targetpts[which.min(targetpts)],1]) & (tuned[1] <= target[targetpts[which.max(targetpts)],1])) 
     }  
          
# final summary plots
 if(genplot)
  {
    par(mfrow=c(3,1))
    if (plotype == 1 || plotype == 2) 
         {
           plot(dat, main="Final stratigraphic correlation points",xlim=c(xmin,xmax),bty="n",lwd=2,cex.axis=1.3,cex.lab=1.3,cex.main=1.4,cex=0.5)
           if(plotype == 1) lines(dat,col="red")
           points(dat[datpts,1],dat[datpts,2],pch=19,col="blue",cex=1.5*size)
           text(dat[datpts,1],dat[datpts,2],labels=1:itick,cex=0.8*size,col="white",font=2)
           plot(target, main="Final astronomical target correlation points",xlim=c(tmin,tmax),bty="n",lwd=2,cex.axis=1.3,cex.lab=1.3,cex.main=1.4,cex=0.5)
           if(plotype == 1) lines(target,col="red")
           points(target[targetpts,1],target[targetpts,2],pch=19,col="blue",cex=1.5*size)
           text(target[targetpts,1],target[targetpts,2],labels=1:ijtick,cex=0.8*size,col="white",font=2)
           plot(tuned, main="Tuned stratigraphic series",bty="n",lwd=2,cex.axis=1.3,cex.lab=1.3,cex.main=1.4,cex=0.5)
           if(plotype == 1) lines(tuned,col="red")
           points(tuned[(datpts-remov),1],tuned[(datpts-remov),2],pch=19,col="blue",cex=1.5*size)
           text(tuned[(datpts-remov),1],tuned[(datpts-remov),2],labels=1:itick,cex=0.8*size,col="white",font=2)
         }
    if (plotype == 3) 
         {
           plot(dat, type="l", main="Final stratigraphic correlation points",xlim=c(xmin,xmax),bty="n",lwd=2,cex.axis=1.3,cex.lab=1.3,cex.main=1.4,cex=0.5)
           points(dat[datpts,1],dat[datpts,2],pch=19,col="blue",cex=1.5*size)
           text(dat[datpts,1],dat[datpts,2],labels=1:itick,cex=0.8*size,col="white",font=2)
           plot(target, type="l", main="Final astronomical target correlation points",xlim=c(tmin,tmax),bty="n",lwd=2,cex.axis=1.3,cex.lab=1.3,cex.main=1.4,cex=0.5)
           points(target[targetpts,1],target[targetpts,2],pch=19,col="blue",cex=1.5*size)
           text(target[targetpts,1],target[targetpts,2],labels=1:ijtick,cex=0.8*size,col="white",font=2)
           plot(tuned, type="l", main="Tuned stratigraphic series",bty="n",lwd=2,cex.axis=1.3,cex.lab=1.3,cex.main=1.4,cex=0.5)
           points(tuned[(datpts-remov),1],tuned[(datpts-remov),2],pch=19,col="blue",cex=1.5*size)
           text(tuned[(datpts-remov),1],tuned[(datpts-remov),2],labels=1:itick,cex=0.8*size,col="white",font=2)
         }

# generate additional summary plots (time-depth map and sedimentation rates
    dev.new(title=paste("Time-space Map and Sedimentation Rates"),height=8,width=4)
    par(mfrow=c(2,1))
    plot(controlPts,bty="n",main="Time-Space Map"); lines(controlPts,col="red")
# calculate sedimentation rates
    numC = length(targetpts)
    thick = controlPts[2:numC,1]-controlPts[1:(numC-1),1]
    dur = controlPts[2:numC,2]-controlPts[1:(numC-1),2]
    sedrate = thick/dur
# make sedimentation rate plot. slightly jitter control depths/heights, for plotting purposes only.
    dsed=double( (numC*2) - 2)
    ssed=double( (numC*2) - 2)
# first datum
    dsed[1] = controlPts[1,1]
    ssed[1] = sedrate[1]
    ii = 2
# middle of section  
    if(numC > 2)
     {
       for (i in seq(2, (numC*2) - 3, by=2))
        {
          dsed[i] = controlPts[ii,1] -  controlPts[ii,1]*0.000001
          dsed[i+1] = controlPts[ii,1] +  controlPts[ii,1]*0.000001
          ssed[i] = sedrate[ii-1]
          ssed[i+1] = sedrate[ii]
          ii = ii + 1
        }
     }   
# last datum
    dsed[(numC*2) - 2] = controlPts[numC,1]
    ssed[(numC*2) - 2] = sedrate[numC-1]
    plot(dsed,ssed, type="l",bty='n',ylab="Sedimentation Rate",xlab="Depth/Height",main="Sedimentation Rates",cex.axis=1.1,cex.lab=1.1)

# end genplot section
  }

    cat("\n * Tuning complete.\n")

    if (output == 1) return(tuned)
    if (output == 2) return(controlPts)
    if (output == 3) return(list(controlPts,tuned))

### END function linage
}