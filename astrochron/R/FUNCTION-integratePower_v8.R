### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### integratePower: determine the variance within a given bandwidth, and 
###                 ratio it to the total power in the spectrum
###                 (SRM: December 13, 2013; December 14, 2013; December 17, 2013;
###                   December 18, 2013; December 24, 2013; January 15-16, 2015;
###                   March 6, 2015; March 19, 2015)
###
###########################################################################

integratePower <- function (spec,flow=NULL,fhigh=NULL,fmax=NULL,unity=F,f0=T,xmin=NULL,xmax=NULL,ymin=NULL,ymax=NULL,npts=NULL,pad=NULL,ydir=1,ncolors=100,h=6,w=9,ln=F,genplot=T,verbose=T)
{

if(verbose) cat("\n----- INTEGRATING POWER SPECTRUM -----\n")


 if(unity && verbose)  cat("\n * NOTE: Spectrum will be normalized such that total variance (up to fmax) is unity.\n")

 if(!unity && is.null(npts))
   { 
       cat("\n**** Please specify the number of points (npts) in the processed time series window.\n")
       stop("**** TERMINATING NOW!")
   }
   
 if(!unity && is.null(pad))
   { 
       cat("\n**** Please specify the total padding length used (pad).\n")
       stop("**** TERMINATING NOW!")
   }   
   
   if(!unity && npts != pad) fac=0.5*(pad/npts)
   if(!unity && npts == pad) fac=0.5 
   if(unity) fac=0.5 

# ensure we have a data frame
  spec=data.frame(spec)

# assign frequencies from first column of spec
  freq=spec[,1]

cols=length(spec)
if(cols > 2)
 {
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
  if(is.null(ymin)) ymin = min(loc)
  if(is.null(ymax)) ymax = max(loc)
}

  if(cols == 2) 
  {
    numrec=1
    sp=spec[,2]
  }  

  spPlot=sp
  if(ln) spPlot=log(spPlot)

  numfreq=length(freq)
  if(verbose) cat("\n * Number of spectra to analyze =",numrec,"\n")
  if(verbose) cat(" * Number of frequencies per spectrum =",numfreq,"\n")

  if(is.null(xmin)) xmin = min(freq)
  if(is.null(xmax)) xmax = max(freq)
  if(is.null(fmax)) fmax = max(freq)


  if(genplot) plotit=T
  if(!genplot)
   {
     if(is.null(flow) || is.null(fhigh)) 
       {
         plotit=T
         cat("\n **** NOTE: flow and/or fhigh have not been specified. Interactive plotting will be activated.\n")
       }  
     if(!is.null(flow) && !is.null(fhigh)) plotit=F
   }
   
   if(plotit && numrec>1)
    {   
# spectrogram plot
     dev.new(height=h,width=w)
     mat <- matrix(c(1, 2, 3, 4), nrow = 1, ncol = 4)
     layout(mat, widths = c(3.5, 1.5, 1.5, 1.5))

# use fields library for access to 'tim.colors'
# useRaster=T results in a faster plotting time.
     xlimset=c(xmin,xmax)
     if (ydir == -1) 
       {
# in this case, reset ylim range.
# note that useRaster=T is not a viable option, as it will plot the results backwards, even though the
#  y-axis scale has been reversed!  This option will result in a slower plotting time.
         ylimset=c(ymax,ymin)
         image(freq,loc,spPlot,xlim=xlimset,ylim=ylimset,col = tim.colors(ncolors),xlab="Frequency",ylab="Location",main="")
         if(is.null(flow) || is.null(fhigh)) mtext("Select frequency band by clicking twice on plot",side=3,line=1)
       }
     if (ydir == 1) 
       {
# useRaster=T results in a faster plotting time.
        ylimset=c(ymin,ymax)
        image(freq,loc,spPlot,xlim=xlimset,ylim=ylimset,col = tim.colors(ncolors),useRaster=T,xlab="Frequency",ylab="Location",main="")       
        if(is.null(flow) || is.null(fhigh)) mtext("Select frequency band by clicking twice on plot",side=3,line=1)
       }

     if(is.null(flow) || is.null(fhigh)) 
# interactive plot selection
      {
       cat("\n   PLEASE CLICK ON PLOT TO SELECT TWO FREQUENCIES TO DEFINE POWER INTEGRATION BAND\n")
# Now overlay x-y plot for graphical interface
       par(new=T)
       plot(-1,-1,xlim=xlimset,ylim=ylimset,xaxs="i",yaxs="i",yaxt='n',bty='n',ylab="",xlab="")
       ff = locator(n = 2, type = "p" ,col="white",lwd=2,pch=3)
       flow=min(ff$x)
       fhigh=max(ff$x)
      }
     abline(v=c(flow,fhigh),col="white",lty=3,lwd=2)
   }


   if(plotit && numrec==1)
    {   
# spectrogram plot
     dev.new(height=5,width=7)
     par(mfrow=c(1,1))
     plot(freq,spPlot,type="l",xlim=c(xmin,xmax),xlab="Frequency",ylab="",bty = "n", lwd = 2, cex.axis = 1.1, cex.lab = 1.1)
     if(is.null(flow) || is.null(fhigh)) 
# interactive plot selection
      {
       mtext("Select frequency band by clicking twice on plot",side=3,line=1)
       cat("\n   PLEASE CLICK ON PLOT TO SELECT TWO FREQUENCIES TO DEFINE POWER INTEGRATION BAND\n")
       ff = locator(n = 2, type = "p" ,col="red",lwd=1,pch=3)
       flow=min(ff$x)
       fhigh=max(ff$x)
      }
     abline(v=c(flow,fhigh),col="red",lty=3,lwd=2)
   }


# isolate spectrum up to fmax
   ifreqMax= which( (freq <= fmax) )
# isolate the frequency band that you want to integrate
   ifreq= which( (freq >= flow) & (freq <= fhigh) )

# only count f0 once (use 0.5 since fac will double others to get power from neg freqs)
   if(f0) sp[1]=sp[1]*0.5

# perform integration
  bandPwr=double(numrec)
  totalPwr=double(numrec)
  for (i in 1:numrec)
    {
      if(numrec>1) 
       {
         if(verbose) cat(" * Integrating window=",i,"; Location=",loc[i],"\n")
         sp1=sp[ifreqMax,i]
         sp2 = sp[ifreq,i]
       }
      if(numrec==1)
       {
         sp1=sp[ifreqMax]
         sp2 = sp[ifreq]
        }
# integrate up to fmax
        totalPwr[i]=sum(sp1)/fac
# integrate band of interest
        bandPwr[i]=sum(sp2)/fac
        if(unity) bandPwr[i]=bandPwr[i]/totalPwr[i]
        if(unity) totalPwr[i] = 1
    }


  if(verbose)
   {
    cat("\n   Frequency Band selected for Power Integration=",flow,",",fhigh,"\n")
    cat("\n   Equivalent periods=",1/flow,",",1/fhigh,"\n\n")
   }

  if(numrec > 1)
   {
    out <- data.frame(cbind(loc,bandPwr,totalPwr,bandPwr/totalPwr))
    colnames(out)<-c("Location","Band Power","Total Power","Band/Total Power")

    if(plotit)
     {
      plot(bandPwr,loc,type="l",ylim=ylimset,main="Band Power",ylab="Location",xlab="Cumulative Power")
      plot(totalPwr,loc,type="l",ylim=ylimset,main="Total Power",ylab="Location",xlab="Cumulative Power")
      plot(bandPwr/totalPwr,loc,type="l",ylim=ylimset,main="Band/Total Power",ylab="Location",xlab="Fraction of Power")
     }
    }
    
   if(numrec==1)
    {
      out <- data.frame(cbind(bandPwr,totalPwr,bandPwr/totalPwr))
      colnames(out)<-c("Band Power","Total Power","Band/Total Power")
    }  

    return(out)

### END function integratePower
}
