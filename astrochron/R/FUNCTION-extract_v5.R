### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### extract: extract record from eha or eAsm output (SRM: July 1, 2013; Aug. 8, 2013
###                                      Feb. 13, 2014; Jan. 16, 2015; March 6, 2015)
###                                                                 
###########################################################################

### note that ydir=1 will not work for eAsm, since raster does not allow uneven sample grid- Feb. 13, 2014

extract <- function (spec,get=NULL,xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,h=6,w=4,ydir=1,pl=0,ncolors=100,genplot=T,verbose=T)
{
  
  if(verbose) cat("\n---- EXTRACTING RECORD ----\n")

# ensure we have a data frame
  spec=data.frame(spec)

# assign frequencies from first column of spec
  freq=spec[,1]

  cols=length(spec)
# assign locations for each spectrum (column headers)
  loc=suppressWarnings(as.numeric(substr(colnames(spec[2:cols]),start=2,stop=100)))
# for negative depth/height/time values, "-" has been changed to "."
# this will create NAs. implement modification of fix recommended by Mathieu Martinez
  neg=grepl(".",substr(colnames(spec[2:cols]), start=2,stop=2),fixed=T)
  fixloc=which(neg)
  if(any(neg)) {loc[fixloc]=-1*as.numeric(substr(colnames(spec[(fixloc+1)]),start=3,stop=100))}
# assign specta (amplitude, power, or probability)
  sp=as.matrix( spec[2:cols] )

  numrec=length(loc)
  numfreq=length(freq)
  if(verbose) cat("\n * Number of spectra/windows =",numrec,"\n")
  if(verbose) cat(" * Number of frequencies per spectrum/sedimentation rates =",numfreq,"\n")

# if value specified for 'get', find closest record
  if(!is.null(get)) 
   {
     pts=which.min(abs(loc-get))
     if(verbose) cat("\n * Extracting nearest record =",loc[pts],"\n")
   }

# if record is not specified, use interactive graphical selection
  if(is.null(get))
   {
# use fields library for access to 'tim.colors'

# for plotting
      if(is.null(xmin)) xmin = min(freq)
      if(is.null(xmax)) xmax = max(freq)
      if(is.null(ymin)) ymin = min(loc)
      if(is.null(ymax)) ymax = max(loc)
      if(pl==0) sp_plot=sp
      if(pl==1) sp_plot=log(sp)
      if(pl==2) 
       {
           sp_plot=t(sp)/(apply(sp,2,max))
           sp_plot=t(sp_plot)
       }    
       
# set up device
      dev.new(height=h,width=w)
      par(mfrow=c(1,1))
      xlimset=c(xmin,xmax)

      if (ydir == -1) 
        {
# in this case, reset ylim range.
# note that useRaster=T is not a viable option, as it will plot the results backwards, even though the
#  y-axis scale has been reversed!  This option will result in a slower plotting time.
           ylimset=c(ymax,ymin)
           image(freq,loc,sp_plot,xlim=xlimset,ylim=ylimset,col = tim.colors(ncolors),xlab="Frequency",ylab="Location",main="")
           mtext("Select record to extract",side=3,line=2)
           mtext("by clicking on symbol at right of plot",side=3,line=1)           
         }

      if (ydir == 1) 
        {
# useRaster=T results in a faster plotting time.
           ylimset=c(ymin,ymax)
           image(freq,loc,sp_plot,xlim=xlimset,ylim=ylimset,col = tim.colors(ncolors),useRaster=T,xlab="Frequency",ylab="Location",main="")       
           mtext("Select record to extract",side=3,line=2)
           mtext("by clicking on symbol at right of plot",side=3,line=1)
        }

## this script modified from '?identify' in R
identifyPch <- function(x, y=NULL, n=length(x), pch=19, ...)
{
    xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
    sel <- rep(FALSE, length(x)); res <- integer(0)
    while(sum(sel) < n) {
        ans <- identify(x[!sel], y[!sel], n=1, plot=F, ...)
        if(!length(ans)) break
        ans <- which(!sel)[ans]
        points(x[ans], y[ans], pch = pch, col="black")
        sel[ans] <- TRUE
        res <- c(res, ans)
    }
    res
}

      par(new=T)
#      putx=(xmin+xmax)/2
      putx=xmax+(xmax*.02)
      plot(rep(putx,length(loc)),loc,xlim=xlimset,ylim=ylimset,xaxs="i",yaxs="i",yaxt='n',bty='n',ylab="",xlab="",pch=1,cex=.5,col="black",xpd = TRUE)
      cat("\n ---- Interactively select record ----\n")
      cat("\n *****  Select by clicking on point  *****\n")
      pts <- identifyPch(rep(putx,length(loc)),loc,n=1)
      abline(h=loc[pts],col="white",lwd=2,lty=4)
      if(verbose) cat("\n * Extracting record=",loc[pts],"\n")
# end is.null(get) section
   } 

   out=data.frame(cbind(freq,sp[,pts]))

# Now plot
  if(genplot)
    {
      dev.new(height=4,width=7)
      par(mfrow=c(1,1))
      if(pl==0) plot(out,type="l",col="red",xlim=c(xmin,xmax),cex.axis=1.1,cex.lab=1.1,xlab="",ylab=" ",main="Extracted Result",bty="n",lwd=2)
      if(pl==1) plot(out,type="l",col="red",log="y",xlim=c(xmin,xmax),cex.axis=1.1,cex.lab=1.1,xlab="",ylab=" ",main="Extracted Result",bty="n",lwd=2)
      if(pl==2) plot(out[,1],out[,2]/max(out[,2]),type="l",col="red",xlim=c(xmin,xmax),cex.axis=1.1,cex.lab=1.1,xlab="",ylab=" ",main="Extracted Result",bty="n",lwd=2)
# end genplot section
    }

  return(out)

### END function extract
}