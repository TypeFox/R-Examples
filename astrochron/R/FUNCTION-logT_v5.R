### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### function logT - (SRM: January 30, 2012; December 12, 2012; April 24, 2013
###                       May 20, 2013; June 16, 2015; October 8, 2015)
###
### apply log transformation
###########################################################################

logT <- function (dat,c=0,opt=1,genplot=T,verbose=T) 
{

   if(verbose) cat("\n----- PERFORMING log TRANSFORM OF STRATIGRAPHIC SERIES-----\n")
   dat <- data.frame(dat)
# a is a constant to add before log transformation
   npts <- length(dat[,1]) 
   if(verbose) cat(" * Number of data points=", npts,"\n")
   ii=length(dat)
   if(verbose) cat(" * Number of columns=", ii,"\n")
   ij=ii
   if(ii<=2) ij = 1
   if(ii>2)  ij = ii - 1
   trans <- double(npts*ij)
   dim(trans)<-c(npts,ij)
   for (j in 1:ij)
    {
      if(ii==1) x=dat[,1]
      if(ii>=2) x=dat[,j+1]
      if(opt==1) trans[,j] <- log(x + c )
      if(opt==2) trans[,j] <- log10(x + c )
    }  
 
   if(any(!is.finite(trans)))
    {
       if(verbose) cat("***** WARNING: NA, Nan, -Inf or Inf values produced by transformation!\n")
    } 
    
   if(ii==1) out <- data.frame(trans)
   if(ii>=2) out <- data.frame(cbind(dat[,1],trans))

   if(genplot && ii <=2)
     {
### plot data series. Note, cex is the factor by which to increase or decrease default symbol size
      par(mfrow=c(2,2))
      if(ii==1)
        {
          plot(out[,1], cex=.5,xlab="Location",ylab="Transformed Value", main="Transformed Series")
          lines(out[,1])
        }  
      if(ii==2)
        {
          plot(out, cex=.5,xlab="Location",ylab="Transformed Value",main="New Stratigraphic Series")
          lines(out)
        }  
### plot the denisty and the histogram together
      hist(trans,freq=F) 
      lines(density(trans, bw="nrd"),col="red"); grid()
### boxplot
      boxplot(trans)
### Normal probabilty plot (Normal Q-Q Plot)
      qqnorm(trans); qqline(trans, col="red");grid()
     }
   
   return(out)

### END function logT
}
