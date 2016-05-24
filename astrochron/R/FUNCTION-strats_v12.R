### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### strats - (SRM: January 25, 2012; March 9, 2012; April 24, 2013; 
###                May 20, 2013; August 6, 2013; August 9, 2013; 
###                July 31, 2014; August 1, 2014; August 12-14, 2014;
###                November 7, 2014; April 7-8, 2015; May 18, 2015;
###                September 10, 2015; October 17, 2015)
### 
### summary statistics for stratigraphic data series
###########################################################################

strats <- function (dat,output=0,genplot=1)
{

   cat("\n----- DETERMINING SUMMARY STATISTICS FOR STRATIGRAPHIC SERIES -----\n")
   dat <- data.frame(dat)
   npts <- length(dat[,1])    
   cat(" * Number of data points=",npts,"\n")

   ii=length(dat)
   
### set up plots   
   if(ii==1) par(mfrow=c(2,2)) 
   if(ii>=2) par(mfrow=c(2,3))
   
   if(ii>2) 
     {
       cat("\n**** WARNING: Multiple variables have been detected.\n")
       cat("              Will only evaluate first variable (in second column).\n") 
     }  

### check for NAs
   rmNA=F
   if( any(is.na(dat[,1])) ) 
    { 
     rmNA=T
     if(ii==1) dat = data.frame(subset(dat[,1], !is.na(dat[, 1])))
     if(ii>=2) dat = subset(dat, !is.na(dat[, 1]))
    } 
   if( ii>=2 && any(is.na(dat[,2])) )     
    { 
     rmNA=T
     dat = subset(dat, !is.na(dat[, 2]))
    }

   if(rmNA) 
    {
      cat("\n**** WARNING: NA values detected, and will be removed.\n\n")
      npts <- length(dat[,1])  
    }

   cat(" * New number of data points=",npts,"\n")

   if(ii>=2)
   {

### evaluate sampling statistics
   t1<-dat[1:(npts-1),1]
   t2<-dat[2:(npts),1]
   dt=t2-t1
      
### check for sorting   
   srt=F
   if(min(dt)<0 && max(dt)>0) srt=T
   
   if(srt)
    {
      cat("\n**** WARNING: Stratigraphic series is not sorted.\n")
      cat(" * Sorting data series into increasing order.\n\n")
      dat <- dat[order(dat[1], na.last = NA, decreasing = F),]
      npts <- length(dat[,1]) 
### now reevaluate sampling statistics
      t1<-dat[1:(npts-1),1]
      t2<-dat[2:(npts),1]
      dt=t2-t1   
    }
   
   cat(" * Duration of stratigraphic series=",abs(dat[npts,1]-dat[1,1]),"\n")
   cat("     Minimum=",min(dat[,1]),"\n")
   cat("     Maximum=",max(dat[,1]),"\n")

   dt=abs(dt)
   dtMin=min(dt) 
   dtMax=max(dt)
   dtMean=mean(dt)     
   dtMedian=median(dt)
   
   cat("\n * SAMPLING INTERVAL / EQUIVALENT NYQUIST FREQUENCY:\n")
   cat("     Mean=", dtMean," / " ,1/(2*dtMean),"\n")
   cat("     Median=",dtMedian," / ",1/(2*dtMedian),"\n")
   cat("     Maximum=",dtMax," / ",1/(2*dtMax),"\n")
   cat("     Minimum=", dtMin," / ",1/(2*dtMin),"\n")
   cat(" * Variance of sampling interval=",var(dt),"\n")

# generate cumulative plot data for dt
# versus % of data
   dtSort=sort(dt)
   psum=100*(1:length(dt)/length(dt))
# versus % of time/space
   csum=cumsum(dtSort)
   psum2=100*csum/max(csum)

   if(genplot > 0)
    {
      epsm=1e-9      
### do not output summary plots of dt if evenly sampled   
      if( dtMax-dtMin < epsm ) par(mfrow=c(2,2))
### output summary plots of dt if unevenly sampled       
      if( dtMax-dtMin > epsm ) 
       {
### now output summary plots of dt
         plot(t1,dt, cex=.5, xlab="location",ylab="dt",main="dt by location")
         lines(t1,dt)
### cumulative dt plot
         if(genplot==1)
          {
            plot(psum,dtSort,cex=.5, xlab="Percent",ylab="dt",main="cumulative dt plot")
            lines(psum,dtSort)
            points(psum2,dtSort,cex=.5,col="red")
            lines(psum2,dtSort,col="red")
            mtext("black=%points; red=%duration",side=3,cex=0.75,font=2)
          }  
         if(genplot==2)
          {
### plot the denisty and the histogram together
           hist(dt,freq=F, xlab="dt",main="histogram of dt"); lines(density(dt, bw="nrd"),col="red"); grid() 
          }  
       }  
    }
       
### statistics on data values
    cat("\n * Mean data value=",mean(dat[,2]),"\n")
    cat(" * Median data value=",median(dat[,2]),"\n")
    cat(" * Minimum data value=",min(dat[,2]),"\n")
    cat(" * Maximum data value=",max(dat[,2]),"\n")
    cat(" * Variance of data values=",var(dat[,2]),"\n")

    if(genplot>0)
     {
### now plot y variable of data series. Note, cex is the factor by which to increase or decrease default symbol size
      plot(dat[,1],dat[,2], cex=.5, xlab="location",ylab=colnames(dat)[2],main="data")
      lines(dat[,1],dat[,2])
### plot the denisty and the histogram together
      hist(dat[,2],freq=F, xlab=colnames(dat)[2],main="distribution of values"); lines(density(dat[,2], bw="nrd"),col="red"); grid()
### boxplot
      boxplot(dat[,2], ylab=colnames(dat)[2],main="boxplot of values")
### Normal probabilty plot (Normal Q-Q Plot)
      qqnorm(dat[,2]); qqline(dat[,2], col="red"); grid()
    }

### end ii>=2 section
    }

    if(genplot>0 && ii==1)
     {
      plot(dat[,1], cex=.5, xlab="location",ylab=colnames(dat)[1],main="data")
      lines(dat[,1])
### plot the denisty and the histogram together
      hist(dat[,1],freq=F, xlab=colnames(dat)[1],main="distribution of values"); lines(density(dat[,1], bw="nrd"),col="red"); grid()
### boxplot
      boxplot(dat[,1], ylab=colnames(dat)[1],main="boxplot of values")
### Normal probabilty plot (Normal Q-Q Plot)
      qqnorm(dat[,1]); qqline(dat[,1], col="red"); grid()
    }

    if(output==1) return(data.frame(cbind(psum,dtSort)))
    if(output==2) return(data.frame(cbind(psum2,dtSort)))
    if(output==3) return(data.frame(cbind(t1,dt)))
   
### END function strats
}
