### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### function trough : find minima of troughs in series, report those that
###                   fall below a threshold value (February 3, 2015;
###                   February 7, 2015; August 17, 2015)
###
###########################################################################

trough <- function (dat,level=NULL,genplot=T,verbose=T) 
{

if(verbose) cat("\n----- FINDING MINIMA OF TROUGHS, FILTERING AT THRESHOLD VALUE -----\n")

dat <- data.frame(dat)
npts <- as.integer(dim(dat)[1])
ncols <- as.integer(dim(dat)[2])
if(verbose) cat(" * Number of data points=", npts,"\n")
if(verbose) cat(" * Number of columns=", ncols,"\n")

if(ncols == 1) { y <- dat[,1] }
if(ncols == 2) { x <- dat[,1]; y <- dat[,2] }

# FORTRAN wrapper
troughID <- function (npts,y) 
 { 
    F_dat = .Fortran('trough_r',PACKAGE='astrochron',
    
    npts=as.integer(npts),y=as.double(y),
    
    loc=integer(as.integer(npts)),iplat=integer(as.integer(npts)),
    numtrough=integer(1),numplat=integer(1)
    )    
# return the results
    return(F_dat)
 }

if(verbose) cat(" * Identifying minima of troughs\n")
# identify minima of troughs
res <- troughID(npts,y)

numtrough <- res$numtrough
numplat <- res$numplat
loc <- res$loc[1:numtrough]
iplat <- res$iplat[1:numplat]

if(verbose) cat(" * Number of troughs identified=",numtrough,"\n")

if(numplat>0) 
 {
   if(verbose) cat(" * Number of plateau points detected=",numplat,"\n")
   if(ncols == 1)
    {
      plats = data.frame( cbind(iplat,y[iplat]))
      colnames(plats)[1] <- 'Index'
      colnames(plats)[2] <- 'Plateau_Value'
    }
   
   if(ncols == 2)
    {
      plats = data.frame( cbind(iplat,x[iplat],y[iplat]))
      colnames(plats)[1] <- 'Index'
      colnames(plats)[2] <- 'Plateau_Location'
      colnames(plats)[3] <- 'Plateau_Value'
    }
# this warning should be output regardless of whether verbose selected.
      cat("\n**** WARNING: The following plateaus were not evaluated!:\n")
      print(plats)  
 }
  
if(ncols == 1) { ymin <- y[loc] }
if(ncols == 2) { xloc <- x[loc] ; ymin <- y[loc] }

if(numtrough > 0)
 {

if( is.null(level) ) 
 {
    if(verbose) cat("\n * No filtering of troughs applied.\n")
    filt_loc <- loc
    if(ncols == 2) filt_xloc <- xloc
    filt_ymin <- ymin
  }    

skip=F
      
if( !is.null(level) && level < min(ymin))
 {
# this warning should be output regardless of whether verbose selected.
    cat("\n**** WARNING: threshold level for filtering is less than minimum value observed \n")
# in this case, reassign numtrough to zero so nothing is returned
    numtrough=0
    skip=T
 }

if(!skip)
{
      
if( !is.null(level) && level >= min(ymin))
 {
   if(verbose) cat(" * Filtering troughs at threshold of",level,"\n")
   filt_loc <- loc[ymin <= level]
   if(ncols == 2) filt_xloc <- xloc[ymin <= level]
   filt_ymin <- ymin[ymin <= level]
   numtrough=length(filt_ymin)
   if(verbose) cat(" * Number of troughs <=",level,":",numtrough,"\n")
 }
  
if(numtrough>0 && ncols == 1) 
  {
        out = data.frame( cbind(filt_loc,filt_ymin) )
        colnames(out)[1] = 'ID'
        colnames(out)[2] = 'Trough_Value'
  }
if(numtrough>0 && ncols == 2)
  {
        out = data.frame( cbind(filt_loc,filt_xloc,filt_ymin) )
        colnames(out)[1] = 'ID'
        colnames(out)[2] = 'Location'
        colnames(out)[3] = 'Trough_Value'
  }
 
if(numtrough>0 && genplot)
 {
   par(mfrow=c(1,1))
   nfilt=length(filt_loc)
   if(ncols == 1)
     {
      plot(1:npts,y,cex=0.5,main="Data with Trough Minima Identified",xlab="Point Number",ylab="Value",bty="n")
      lines(1:npts,y,col="forestgreen")
      abline(v=filt_loc,col="blue",lty=22)
      points(filt_loc,filt_ymin,pch=1,col='blue')
      if(numtrough >1)
        {
          mtext(filt_loc[seq(1,nfilt,by=2)], side=3,line=0.25,at=filt_loc[seq(1,nfilt,by=2)],cex=0.5,font=4,col="blue")
          mtext(filt_loc[seq(2,nfilt,by=2)], side=3,line=-0.25,at=filt_loc[seq(2,nfilt,by=2)],cex=0.5,font=4,col="blue")
        }
      if(numtrough ==1) mtext(filt_loc, side=3,line=0.25,at=filt_loc,cex=0.5,font=4,col="blue")
    
      if(numplat > 0)
       {
         if(numplat == 1) mtext(iplat, side=3,line=0.25,at=iplat,cex=0.5,font=4,col="red")
         if(numplat > 1) mtext(iplat[seq(1,max(iplat),by=2)], side=3,line=0.25,at=iplat[seq(1,max(iplat),by=2)],cex=0.5,font=4,col="red")
         points(iplat,y[iplat],pch=19,col='red')
       } 
    }   
   if(ncols == 2)
    {
      plot(dat,cex=0.5,main="Data with Trough Minima Identified",xlab="Location",ylab="Value",bty="n")
      lines(dat,col="forestgreen")
      abline(v=filt_xloc,col="blue",lty=22)
      points(filt_xloc,filt_ymin,pch=1,col='blue')
      if(numtrough >1)
        {
          mtext(filt_xloc[seq(1,nfilt,by=2)], side=3,line=0.25,at=filt_xloc[seq(1,nfilt,by=2)],cex=0.5,font=4,col="blue")
          mtext(filt_xloc[seq(2,nfilt,by=2)], side=3,line=-0.25,at=filt_xloc[seq(2,nfilt,by=2)],cex=0.5,font=4,col="blue")
        }
       if(numtrough==1) mtext(filt_xloc, side=3,line=0.25,at=filt_xloc,cex=0.5,font=4,col="blue")
        
      if(numplat > 0)
       {
         if(numplat==1) mtext(iplat, side=3,line=0.25,at=x[iplat],cex=0.5,font=4,col="red")
         if(numplat>1) mtext(iplat[seq(1,length(iplat),by=2)], side=3,line=0.25,at=x[iplat[seq(1,length(iplat),by=2)]],cex=0.5,font=4,col="red")
         points(x[iplat],y[iplat],pch=19,col='red')
       } 
    }   
# end genplot section
 }
}

# end skip=F
}
 
if(numtrough>0) return( out )

### END function trough
}

