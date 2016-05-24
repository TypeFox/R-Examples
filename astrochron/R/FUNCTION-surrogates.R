### This function is a component of astrochron: An R Package for Astrochronology
### Copyright (C) 2015 Stephen R. Meyers
###
###########################################################################
### function ebiSurrogate : generate random phase surrogates as in Ebisuzaki W. (1997),
###       A method to estimate the statistical significance of a correlation 
###       when the data are serially correlated. J Climate, 10, 2147-2153.
###                                 (SRM: Dec 6-11, 2012; May 28-29, 2013; 
###                                       June 5, 2013; February 4, 2015)
###           								
### This is an R-translation of the Matlab code by V. Moran, with modifications 
###  and additional features.
### The code by V. Moran can be found here:
###  http://www.mathworks.com/matlabcentral/fileexchange/10881-weaclim/content/ebisuzaki.m
### See also http://www.ftp.cpc.ncep.noaa.gov/wd51we/random_phase/ 
###########################################################################

surrogates <- function (dat,nsim=1,preserveMean=T,std=T,genplot=T,verbose=T)
{

if(verbose) cat("\n----- GENERATING PHASE-RANDOMIZED SURROGATES -----\n")

dat<- data.frame(dat)
n=length(dat[,1])
if (verbose) cat(" * Number of data points=", n,"\n")
if (verbose) cat(" * Number of simulations=", nsim,"\n")

n2=floor(n/2)
# allow input of data series with one or two columns (if two, second is used)
if(length(dat)==1) x=dat[,1]
if(length(dat)==2) x=dat[,2]
# remove mean, store
xmean=mean(x)
x=x-xmean
  
# calculate fft
y=fft(x)
# amplitude estimate
mod=abs(y)

# initialize matrix for simulation results
# if even number of data points
if (n/2==n2) 
  {
    X <- double(n*nsim)
    dim(X) <- c(n,nsim)
  } else
# if odd number of data points no f(Nyq)
  {
    X <- double((2*n2+1)*nsim)
    dim(X) <- c((2*n2+1),nsim)  
  } 

# start simulation loop
for (i in 1:nsim)
 {

# for even number of data points
   if (n/2==n2) 
     {
# phase randomization using uniform distribution [0,2pi) 
# note from runif: "runif will not generate either of the extreme values unless 
#                    max = min or max-min is small compared to min, and in particular 
#                    not for the default arguments."
#   (f(0) not required, since mean removed)
        theta=runif(n2-1,min=0,max=2*pi)
# assign f(0), which is real
	    ph=0
# assign up to (but not including) f(Nyq)
	    ph=append(ph,theta)	
# determine f(Nyq), which is real
        thetaNyq=sqrt(2)*cos(runif(1,min=0,max=2*pi))
	    ph=append(ph,thetaNyq)
# assign negative frequecies in decreasing order
#  phase should be opposite sign
	    ph=append(ph,rev(-1*theta))     
        recf=mod*exp(1i*ph)
#  adjust f(Nyq), as in Ebisuzaki (1997)
#        recf[n2+1] = y[n2+1]
        recf[n2+1] = mod[n2+1]*thetaNyq
# inverse fft
        X[,i]=Re(fft(recf,inverse=T))/n
        if(std) X[,i]=X[,i]*sd(x)/sd(X[,i])
     } else 
    
# if there are an odd number of frequencies, f(Nyq) does not exist        
     {
# phase randomization using uniform distribution [0,2pi)
        theta=runif(n2,min=0,max=2*pi)
# assign f(0), which is real
	    ph=0
# assign up to (but not including) f(Nyq)
	    ph=append(ph,theta)	
# assign negative frequecies in decreasing order
#  phase should be opposite sign
	    ph=append(ph,rev(-1*theta))       
        recf=mod*exp(1i*ph)
# inverse fft
        X[,i]=Re(fft(recf,inverse=T))/n
        if(std) X[,i]=X[,i]*sd(x)/sd(X[,i])
      }
 }
    
if(preserveMean) X = X + xmean

if(nsim == 1) 
  {
    if(genplot)
     {
### plots
       par(mfrow=c(2,2))
       if(length(dat)==2) 
         {
           plot(dat[,1],X,cex=0.5,xlab="Location",ylab="Surrogate Value",main="Phase-randomized Surrogates")
           lines(dat[,1],X)
          } 
      if(length(dat)==1) 
         {
           plot(1:n,X,cex=0.5,xlab="Location",ylab="Surrogate Value",main="Phase-randomized Surrogates")
           lines(1:n,X)
          }    
### plot the denisty and the histogram together
       hist(X,freq=F,xlab="Surrogate Value",main="Distribution of Surrogates"); lines(density(X, bw="nrd"),col="red"); grid()
### boxplot
       boxplot(X,ylab="Surrogate Value",main="Boxplot of Surrogates")
### Normal probabilty plot (Normal Q-Q Plot)
       qqnorm(X); qqline(X, col="red")
     } 
 
    if(length(dat)==2) return(data.frame( cbind(dat[,1],X) ) )
    if(length(dat)==1) return(data.frame(X))
  }

if(nsim > 1) return(X) 

}
# end function surrogates