#########################################################################
#       $Log: Util.S,v $
#
#               (c) Copyright  1997
#                          by                                   
#      Author: Rene Carmona,Bruno Torresani,Wen-Liang Hwang,Andra Wang   
#                  Princeton University
#                  All right reserved                           
#########################################################################





smoothwt <- function(modulus, subrate, flag = FALSE)
#########################################################################
#      smoothwt:   
#      --------
#       smooth the wavelet (or Gabor) transform in the b direction
#
#       input:
#       ------
#	 modulus: modulus of the wavelet (or Gabor) transform
#	 subrate: size of the smoothing window
#	 flag: if set to TRUE, subsamples the smoothed transform by
#		subrate.	
#
#       output:
#       -------
#        smoothed transform (2D array)
#
#########################################################################
{
  sigsize <- dim(modulus)[1]
  nscale <- dim(modulus)[2]

  if(flag) {
    smodulus <- matrix(0,sigsize,nscale)
    dim(smodulus) <- c(sigsize * nscale, 1)
  }
  else {
    modsize <- as.integer(sigsize/subrate) + 1
    smodulus <- matrix(0,modsize,nscale)
    dim(smodulus) <- c(modsize * nscale, 1)
  }
	
  dim(modulus) <- c(sigsize * nscale, 1)

  z <- .C("Ssmoothwt",
          sm = as.double(smodulus),
          as.double(modulus),
          as.integer(sigsize),
          as.integer(nscale),
          as.integer(subrate),
          as.integer(flag),
          PACKAGE="Rwave")

  if(flag)
   dim(z$sm) <- c(sigsize, nscale)
  else 	
   dim(z$sm) <- c(modsize, nscale)
  z$sm
}

smoothts <- function(ts, windowsize)
#########################################################################
# smooth a time series by averaging window
#
# Arguments: 
#	<ts> time series
#       <windowsize> window size
#########################################################################
{
   sigsize <- length(ts)
   sts <- numeric(sigsize)
   nscale <- 1   
   stssize <- sigsize
   dim(sts) <- c(sigsize, 1)
   dim(ts) <- c(sigsize, 1)
   subrate <- windowsize # subrate <- 1
   
   z <- .C("Smodulus_smoothing",
           as.double(ts),
           sts = as.double(sts),
           as.integer(sigsize),
           stssize = as.integer(stssize),            
           as.integer(nscale),
           as.integer(subrate),
           PACKAGE="Rwave")

  sts <- z$sts
  sts
}

adjust.length <- function(inputdata)
#########################################################################
# adjust.length adds zeros to the end of the data if necessary so that 
# its length is a power of 2.  It returns the data with zeros added if
# nessary and the length of the adjusted data.
#
# Arguments: 
#	<inputdata> is either a text file or an S object containing
# 		    data.
#########################################################################
{
  if (is.character(inputdata))
    s <- scan(inputdata)
  else
    s <- inputdata
  np <- length(s)
  
  pow <- 1
  while ( 2*pow < np )
    pow <- 2*pow
  new.np <- 2*pow
  
  if ( np == new.np )
    list( signal=s, length=np )
  else {
    new.s <- 1:new.np
    new.s[1:new.np] <- 0
    new.s[1:np] <- s
    list( signal=new.s, length=new.np )
  }
}

#########################################################################
#      npl:
#      ---
#       plots with n rows.
#
#       input:
#       ------
#	 nrow: number of rows
#
#########################################################################
npl <- function(nbrow)
{
  par(mfrow = c(nbrow, 1))
  cat("")
}

#########################################################################
#      SpecGen:
#      --------
#       Estimate power spectrum
#
#       input:
#       ------
#	 input: input signal.
#        spans: length of Daniell's smoother.
#        taper: relative size of tapering window (cosine taper).
#
#########################################################################
SpecGen <- function(input, spans = 9, taper = 0.2)
{
  tmp <- spec.pgram(input,spans,taper)
  tmp1 <- tmp$spec
  tmp <- 10^(tmp1/10)
  spec <- tmp
  spec
}

#########################################################################
#      SampleGen:
#      ----------
#       generate a real sample with prescribed spectrum      
#
#       input:
#       ------
#	 nspec: spectral density
#
#########################################################################
SampleGen <- function(spec)
{
  sqspec <- sqrt(spec)
  d <- length(sqspec)
  size <- d 
  tmp <- rnorm(2 * size)
  real <- numeric(size)
  imag <- numeric(size)
  real <- tmp[1:size]
  a <- size+1
  b <- 2 * size
  imag <- tmp[a:b]
  real1 <- real * sqspec[1:size]
  imag1 <- imag * sqspec[1:size]
  real2 <- numeric(2 * size)
  imag2 <- numeric(2 * size)
  real2[1:size] <- real1
  imag2[1:size] <- imag1
  for(i in 2:size) {
    real2[b-i+2] <- real1[i]
    imag2[b-i+2] <- -imag1[i]
  }	
  imag2[1] <- 0
  tmp1 <-  fft(real2 + 1i*imag2, inverse=TRUE)
  sample <- Re(tmp1)/sqrt((2*size))
  sample
}

hurst.est <- function(wspec, range, nvoice, plot=TRUE)
#########################################################################
#      hurst.est:   
#      ----------
#       estimate Hurst exponent from (part of) wavelet spectrum
#
#       input:
#       ------
#	 wspec: wavelet spectrum.
#	 nvoice: number of voices.
#	 plot: if set to TRUE, displays regression line on current plot.
#
#       output:
#       -------
#        intercept and slope of regression line. Slope is Hurst exponent.
#
#########################################################################
{
  loctmp <- lsfit(range,log(wspec[range],base=2^(2/nvoice)))
  if(plot) abline(loctmp)
  loctmp$coef
}

wspec.pl <- function(wspec, nvoice)
#########################################################################
#      wspec.pl
#      ----------
#       Displays normalized log of wavelet spectrum
#
#       input:
#       ------
#	 wspec: wavelet spectrum.
#	 nvoice: number of voices.
#
#########################################################################
{
  plot.ts(log(wspec, base=2^(2/nvoice)))
  title("log(wavelet spectrum)", xlab="log(scale)", ylab="V(a)")
  cat(" ")
}
