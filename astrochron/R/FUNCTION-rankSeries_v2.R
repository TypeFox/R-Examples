### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### rankSeries: generate rank stratigraphic series from bedding thickness
###              data, then interpolate to specified interval.
###              (June 25, 2014; July 21, 2015)
###
###########################################################################

rankSeries <- function (dat,dt=NULL,genplot=T,verbose=T)
{

  if (verbose) cat("\n----- CREATE LITHOFACIES RANK SERIES -----\n")
  dat<-data.frame(dat)
  ipts <- length(dat[,1]) 
  if (verbose) cat(" * Number of beds input=", ipts,"\n")

# add up bed thicknesses
  height=append(0,cumsum(dat[,1]))
                 
# add and subract small amount from each datum
      small=0.0000001
      heightLo=double(ipts)
      heightUp=double(ipts)
      heightLo=height[2:(ipts+1)]-small
      heightUp=height[2:(ipts+1)]+small

# first datum      
      heightNew=double((2*ipts))
      rank=double((2*ipts))
      heightNew[1]=height[1]
      rank[1]=dat[1,2]
      j=2
      for(i in 2:(ipts))
       {
          heightNew[j]=heightLo[i-1]
          rank[j]=dat[i-1,2]
          heightNew[j+1]=heightUp[i-1]
          rank[j+1]=dat[i,2]
          j=j+2
       }
# last datum    
       heightNew[j]=height[ipts+1]
       rank[j]=dat[ipts,2]

      if(is.null(dt)) 
       {
         dt = min(dat[,1])/5
         if (verbose) cat(" * Interpolating to grid of", dt,"\n")
       }  
      dat2=data.frame(cbind(heightNew,rank))
   
# interpolate to specified sampling interval
      out=linterp(dat2,dt=dt,genplot=F,verbose=F)

if(genplot)
  {
### plots
    par(mfrow=c(2,2))
    plot(out,cex=0.5,xlab="Stratigraphic Height",ylab="Bed Rank",main="Rank series (blue circle=base of bed)", type="l"); points(height[1:(ipts)],dat[,2], col="blue",cex=1); points(out, col="red",cex=0.5)
### plot the denisty and the histogram together
#    hist(out[,2],freq=F,xlab="Bed Rank",main="Distribution of bed rank"); lines(density(out[,2], bw="nrd"),col="red"); grid()
# remove density as it sometimes fails
    hist(out[,2],freq=F,xlab="Bed Rank",main="Distribution of bed rank"); grid()
### boxplot
    boxplot(out[,2],main="Box plot of bed rank")
### Normal probabilty plot (Normal Q-Q Plot)
    qqnorm(out[,2]); qqline(out[,2], col="red")
  }
  
  return(out)

### END function rankSeries
}
