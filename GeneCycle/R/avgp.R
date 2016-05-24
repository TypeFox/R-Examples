### avgp.R  (2004-12-11)
###
###    Average periodogram and related stuff
###
### Copyright 2003-04 Konstantinos Fokianos and Korbinian Strimmer
###
### This file is part of the `GeneCycle' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA



# periodogram  of a single time series x
periodogram.spec.single <- function(x, method="builtin")
{
    if (method=="builtin")
    {
      # demean but do not detrend to avoid artefacts around zero
      spec <- spectrum(x, taper=0, plot=FALSE, fast=FALSE, detrend=FALSE, demean=TRUE)$spec
    } 
    else if (method =="clone")
    {
       x <- x - mean(x) # demean
       
       N<-length(x)
       xfft <- fft(x)
       pgram <- (Mod(xfft)^2/(N))  # squared modulus
       #pgram <- Re(xfft*Conj(xfft)/(N )) #this is the same 
       
       spec <- pgram[1 + (1:floor(N/2)) ] 
    }
    else if (method =="smooth")
     {
      # demean but do not detrend to avoid artefacts around zero
      spec <- spectrum(x, taper=0, plot=FALSE, fast=FALSE, detrend=FALSE, demean=TRUE, span=3)$spec
    } 
   
    return(spec)
}


# periodogram  of multiple time series x
periodogram.spec <- function(x, method="builtin")
{
    f <- periodogram.freq(x)
     
    xm <- as.matrix(x)
     
    num.series <- dim(xm)[2] # number of columns
    spec.matrix <- matrix(NA, nrow=length(f), ncol=num.series)
    for (i in 1:num.series)
    {
       spec.matrix[,i] <- periodogram.spec.single(xm[,i], method=method)
    }
    return(spec.matrix)
}

 
# corresponding frequencies (ranging from 0 to 1/frequency(x))
periodogram.freq <- function(x,  method="builtin")
{
    z <- as.matrix(x)[,1] # use first time series (in first column)
    
    if (method=="builtin")
    {
      # demean but do not detrend to avoid artefacts around zero
      freq <- spectrum(z, taper=0, plot=FALSE, fast=FALSE, detrend=FALSE, demean=TRUE)$freq
    } 
    else if (method =="clone")
    {
       xfreq <- frequency(z)
       N <-length(z)
       Nspec <- floor(N/2)
       freq <- seq(from = xfreq/N, by = xfreq/N, length = Nspec)  
    }
    else if (method =="smooth")
     {
      # demean but do not detrend to avoid artefacts around zero
      freq <- spectrum(z, taper=0, plot=FALSE, fast=FALSE, detrend=FALSE, demean=TRUE, span=3)$freq
    } 
   
    return(freq)
}

# periodogram
periodogram <- function(x,  method="builtin")
{
  list(spec=periodogram.spec(x, method=method), freq=periodogram.freq(x, method=method))
}


# Average Periodogram:
avgp <- function(x, title=deparse(substitute(x)), plot=TRUE, angular = FALSE, ...)
{
    f <- periodogram.freq(x, ...)   
    if (angular) f <- 2*pi*f  # use angular frequencies
    spec.matrix <- periodogram.spec(x, ...)
    avg.spec <- apply(spec.matrix,1,mean)
    out = list(freq=f, avg.spec=avg.spec, title=title)
    
    if (plot)
    {
        plot(out[[1]], out[[2]],  type="l", 
        xlab="Fourier Frequencies", ylab="Average Periodogram", ...)
        title(main=title)
	
        return(invisible(out))
    }
    else return(out)  
}


# returns back the m dominant frequencies in single time series periodogram
dominant.freqs.single <- function(x, m=1, ...)
{    
    spec <- periodogram.spec(x, ...) 
    freq <- periodogram.freq(x, ...)
    
    sorted.freqs <- freq[order(-spec)]
    
    sorted.freqs[1:m]
}


# dito, but now also for multiple time series
dominant.freqs <- function(x, m=1, ...)
{
    xm <- as.matrix(x)
     
    num.series <- dim(xm)[2] # number of columns
    freq.matrix <- matrix(NA, nrow=m, ncol=num.series)
    for (i in 1:num.series)
    {
       freq.matrix[,i] <- dominant.freqs.single(xm[,i], m=m, ...)
    }
    
    return(freq.matrix)
}


