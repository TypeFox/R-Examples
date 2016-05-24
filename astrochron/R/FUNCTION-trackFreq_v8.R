### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### trackFreq: sort EHA output to isolate peaks in specified frequency range 
###             (SRM: June 14, 2013; June 16, 2013; June 26, 2013; July 26, 2013;
###                   December 8, 2013; February 5, 2014; January 16, 2015;
###                   August 17, 2015)
###
###########################################################################

trackFreq <- function (spec,threshold=NULL,pick=T,fmin=NULL,fmax=NULL,dmin=NULL,dmax=NULL,xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,h=6,w=4,ydir=1,ncolors=100,genplot=T,verbose=T)
{
  
  if(verbose) 
    {
      cat("\n----            FREQUENCY DOMAIN MINIMAL TUNING:              ----\n")
      cat("----       TRACK SPATIAL FREQUENCY DRIFT IN EHA PLOTS         ----\n")
      cat("\n----           SORTING EHA OUTPUT TO ISOLATE PEAKS            ----\n")
    }
    
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
  if(verbose) cat("\n * Number of spectra to analyze =",numrec,"\n")
  if(verbose) cat(" * Number of frequencies per spectrum =",numfreq,"\n")

# for analysis
  if(is.null(fmin)) fmin = min(freq)
  if(is.null(fmax)) fmax = max(freq)
  if(is.null(dmin)) dmin = min(loc)
  if(is.null(dmax)) dmax = max(loc)
  
  if(genplot) 
   {
# use fields library for access to 'tim.colors'

# for plotting
      if(is.null(xmin)) xmin = min(freq)
      if(is.null(xmax)) xmax = max(freq)
      if(is.null(ymin)) ymin = min(loc)
      if(is.null(ymax)) ymax = max(loc)

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
           image(freq,loc,sp,xlim=xlimset,ylim=ylimset,col = tim.colors(ncolors),xlab="Frequency",ylab="Location",main="Frequency Domain Minimal Tuning")
         }

      if (ydir == 1) 
        {
# useRaster=T results in a faster plotting time.
           ylimset=c(ymin,ymax)
           image(freq,loc,sp,xlim=xlimset,ylim=ylimset,col = tim.colors(ncolors),useRaster=T,xlab="Frequency",ylab="Location",main="Frequency Domain Minimal Tuning")       
        }
        
# end genplot section
   } 

# isolate portion of the results that you want to analyze (part 1: 'freq')
   ifreq= which( (freq >= fmin) & (freq <= fmax) )
   freq2=freq[ifreq]
   freq2=freq2[freq2 <= fmax]

# perform sorting
  res=rep(NA,numrec)
  for (i in 1:numrec)
    {
      if(verbose) cat(" * Sorting window=",i,"; Location=",loc[i],"\n")
# isolate portion of the results that you want to analyze (part 2: 'sp')
      sp2 = sp[ifreq,i]
      res[i]=peak(cbind(freq2,sp2),level=threshold,genplot=F,verbose=F)[2]
    }


# now rearrange results into two column format (this is potentially slow, and should be optimized!)
  outfreq = 0
  outloc = 0
  for(i in 1:numrec)
   {
# note, is.null sometimes (always?) fails, so use is.na     
     test=(data.frame(res[i]))[1,1]
     if(!is.null(res[i]) || !is.na(test))
      {
        resIn = (data.frame(res[i]))[,1]
        locIn = rep(loc[i],length(resIn))
        outfreq=append(outfreq,resIn)
        outloc=append(outloc,locIn)
      } 
   }    

# remove first 'dummy' value
  outloc=outloc[2:length(outloc)]
  outfreq=outfreq[2:length(outfreq)]
  out=data.frame(cbind(outloc,outfreq))
# Sort to ensure Depth/Height is in increasing order  
#  out <- out[order(out[1],na.last=NA,decreasing=F),]

## this script modified from '?identify' in R
identifyPch <- function(x, y=NULL, n=length(x), pch=19, ...)
{
    xy <- xy.coords(x, y); x <- xy$x; y <- xy$y
    sel <- rep(FALSE, length(x)); res <- integer(0)
    while(sum(sel) < n) {
        ans <- identify(x[!sel], y[!sel], n=1, plot=F, ...)
        if(!length(ans)) break
        ans <- which(!sel)[ans]
        points(x[ans], y[ans], pch = pch, col="white")
        sel[ans] <- TRUE
        res <- c(res, ans)
    }
    res
}

  if(is.null(threshold) && genplot && verbose) 
    {
      cat("\n**** WARNING: By default all peak maxima are reported and plotted.\n") 
      cat("              If there are many peaks, this will cause your plot to be mostly white.\n")
      cat("              Set threshold to cull peaks.\n") 
    }  


# Now plot, identifying location of peaks
  if(genplot)
    {
       par(new=T)
       plot(out[,2],out[,1],xlim=xlimset,ylim=ylimset,xaxs="i",yaxs="i",yaxt='n',bty='n',ylab="",xlab="",pch=1,col="white")
# if interactive point identification selected
       if(pick)
        {
          cat("\n           ---- INTERACTIVELY SELECT FREQUENCIES ----\n")
          cat("\n           *****  Select path by clicking     *****\n")
          cat("   Stop by pressing ESC-key (Mac) or STOP button (Windows)\n")
          pts <- identifyPch(out[,2], out[,1])
          out=data.frame(cbind(out[pts,1],out[pts,2]))
        }  
# end genplot section
    }

       cat("\n * Peaks identified.\n")
       cat("\n * Now:\n")
       cat("   (1) If desired, further cull peaks with function 'idPts' or 'delPts'\n")       
       cat("   (2) Convert the spatial frequencies to sedimentation rates using function 'freq2sedrate'\n")
       cat("   (3) Integrate sedimentation rate curve using function 'sedrate2time'\n")
       cat("   (4) Tune proxy record using function 'tune'\n")

# sort out to ensure increasing depth/height
       out <- out[order(out[1],na.last=NA,decreasing=F),]

       return(out)

### END function trackFreq
}
