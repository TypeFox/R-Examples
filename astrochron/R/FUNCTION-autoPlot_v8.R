### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### autoPlot: automatically plot all data in data frame - VERTICAL plots 
###          (SRM: Oct. 30, 2012; Jan. 18, 2013; May 20, 2013; June 19, 2013;
###                June 28, 2013; August 15, 2013; April 7, 2015)
###
###########################################################################

autoPlot <- function (dat,cols=NULL,nrows=NULL,ydir=-1,smooth=0,xgrid=1,output=F,genplot=T,verbose=T)
{

cat("\n----- PLOTTING (AND SMOOTHING) STRATIGRAPHIC DATA SERIES-----\n")

# ensure we have a data frame
dat <- data.frame(dat)
npts <- length(dat[,1])
if(verbose) cat(" * Number of data points input=", npts,"\n")

# ensure data is sorted to increasing xID order (for ksmooth)
#  missing depths (NA) are removed during sort
if(verbose) cat(" * Sorting into increasing order\n")
dat <- dat[order(dat[1],na.last=NA,decreasing=F),]

smoothScaled= smooth 

if(ydir == 1) ylimset = c(min(dat[,1]), max(dat[,1]))
if(ydir == -1) ylimset = c(max(dat[,1]), min(dat[,1]))

# if cols is not explicitly defined, will use all columns in dat
if(is.null(cols))
  {
    ncols= ( length(dat) - 1 )
    cols = 2:length(dat)
  }
  
if(!is.null(cols))
  {
    ncols=length(cols)
  }
  
if(verbose) cat(" * Number of variables to plot=", ncols,"\n")

if(genplot) 
 {
  par(mar = c(4, 2.5, 1, 2))
  if(is.null(nrows))
   {
    if(ncols<=4) nrows = 1
    if(ncols>4) nrows = ceiling(sqrt(ncols-1))
    ncols1 = nrows
    par(mfrow = c(nrows, ncols1))
   }

  if(!is.null(nrows))
   {
    par(mfrow=c(nrows,ceiling(ncols/nrows)))
   }
 }

# set up smooth
smoothed = rep(NA,npts*ncols)
dim(smoothed) <- c(npts,ncols)  
storename=character(ncols)

for (i in 1:ncols)
 {
  loc=cols[i]
  
  if (smooth == 0) 
    { 
      smoothed[,i] <- dat[,loc]
      storename[i] <- colnames(dat[loc])
      xID <- dat[,1]
     }
  
  if (smooth != 0)
    {
# if x.points = dat[,1], will evalute at original sample locations only
      if(xgrid==1) smooth2 <- ksmooth(dat[,1],dat[,loc],kernel=c("normal"),bandwidth=smoothScaled,x.points=dat[,1])
# if x.points undefined, will evalute on even grid, covering total x range
      if(xgrid==2) smooth2 <- ksmooth(dat[,1],dat[,loc],kernel=c("normal"),bandwidth=smoothScaled)
      smoothed[,i] <- smooth2$y
      storename[i] <- colnames(dat[loc])
      xID <- smooth2$x
     }
     
  colnames(smoothed) <- storename        
  if(genplot) 
   {
     plot(smoothed[,i],xID, cex=0.5, ylim=ylimset,ylab=colnames(dat[1]),xlab=colnames(dat[loc]))
     lines(smoothed[,i],xID, col="black")
   }  
 }

if(output) 
  {
    out <- data.frame(cbind(xID,smoothed))
    colnames(out)[1] <- colnames(dat[1])  
    return(out)
  }
   
### END function autoPlot
}
