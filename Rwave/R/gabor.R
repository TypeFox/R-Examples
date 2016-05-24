#########################################################################
#      $Log: Gabor.S,v $
#########################################################################
#
#               (c) Copyright  1997
#                          by                                   
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  Princeton University 
#                  All right reserved                           
#########################################################################





cgt <- function(input, nvoice, freqstep = (1/nvoice),
                scale = 1, plot = TRUE)
#########################################################################
#       cgt:   
#       ---
#        continuous Gabor transform function:
# 	  compute the continuous Gabor transform with gaussian window.
#
#       input:
#       ------
# 	 input: input signal (possibly complex-valued)
#	 nvoice: number of frequency bands
#        freqstep: sampling rate for the frequency axis
#        scale: size of the window
#	 plot: if set to TRUE, displays the modulus of cwt on the graphic
#		device.
#
#       output:
#       -------
#        output: continuous (complex) gabor transform
#
#########################################################################
{
   oldinput <- input
   isize <- length(oldinput)

   tmp <- adjust.length(oldinput)
   input <- tmp$signal
   newsize <- length(input)

   pp <- nvoice
   Routput <- matrix(0, newsize, pp)
   Ioutput <- matrix(0, newsize, pp)
   output <- matrix(0, newsize, pp)
   dim(Routput) <- c(pp * newsize, 1)
   dim(Ioutput) <- c(pp * newsize, 1)
   dim(input) <- c(newsize, 1)


   z <- .C("Sgabor",
           as.double(input),
           Rtmp = as.double(Routput),
           Itmp = as.double(Ioutput),
           as.integer(nvoice),
           as.double(freqstep),
           as.integer(newsize),
           as.double(scale),
           PACKAGE="Rwave")

   Routput <- z$Rtmp
   Ioutput <- z$Itmp
   dim(Routput) <- c(newsize, pp)
   dim(Ioutput) <- c(newsize, pp)

   output <- Routput[1:isize,] + 1i*Ioutput[1:isize,]
   if(plot) {
      image(1:isize, seq(0, nvoice*freqstep/2, length=nvoice),
            Mod(output), xlab="Time", ylab="Frequency")
      title("Gabor Transform Modulus")
   }
   output
}



vgt <- function(input, frequency, scale, plot = FALSE)
#########################################################################
#      vgt:   
#      ---
#       continuous Gabor transform on one frequency:
# 	 compute the continuous Gabor transform with (complex-valued)
#	 gaussian window
#
#       input:
#       ------
# 	 input: input signal (possibly complex-valued)
#        frequency: value of the frequency
#        scale: size of the window
#	 plot: if set to TRUE, plotss the real part of cgt on the graphic
#		device.
#
#       output:
#       -------
#        Routput + i Ioutput: voice gabor transform (complex 1D array)
#
#########################################################################
{
   oldinput <- input
   isize <- length(oldinput)

   tmp <- adjust.length(oldinput)
   input <- tmp$signal
   newsize <- length(input)

   Routput <- numeric(newsize)
   Ioutput <- numeric(newsize)
   dim(input) <- c(newsize,1)

   z <- .C("Svgabor",
            as.double(input),
            Rtmp = as.double(Routput),
            Itmp = as.double(Ioutput),
            as.double(frequency),
            as.integer(newsize),
            as.double(scale),
           PACKAGE="Rwave")

   Routput <- z$Rtmp
   Ioutput <- z$Itmp
   if(plot==TRUE) {
      plot.ts(Re(z$tmp));
      title("Real part of Gabor transform");
   }

   Routput[1:isize] + 1i*Ioutput[1:isize]
}



gabor <- function(sigsize, location, frequency, scale)
#########################################################################
#       gabor:   
#       ------
#        Generates a Gabor for given location and frequency
#
#       input:
#       ------
#        sigsize: signal size (dimension of the array)
#        location: location of the wavelet
#        frequency: value of the frequency
#        scale: size of the window
#
#       output:
#       -------
#        z$gabor.r + z$gabor.i * i: gabor (complex 1D array
#          of size sigsize)
#
#########################################################################
{
   gabor.r <- numeric(sigsize)
   gabor.i <- numeric(sigsize)


   z <- .C("gabor_time",
            as.double(frequency),
            as.double(scale),
            as.integer(location),
            gabor.r = as.double(gabor.r),
            gabor.i = as.double(gabor.i),
            as.integer(sigsize),
           PACKAGE="Rwave")

   z$gabor.r + 1i*z$gabor.i
}


vecgabor <- function(sigsize, nbnodes, location, frequency, scale)
#########################################################################
#       vecgabor:   
#       --------
#        Generates Gabor functions for given locations and frequencies
#	 on a ridge.
#
#       input:
#       ------
#        sigsize: signal size (dimension of the array)
#        nbnodes: number of ridge samples
#        location: b coordinates of the ridge samples
#        frequency: acoordinates of the ridge samples
#        scale: size of the window
#
#       output:
#       -------
#        z$gabor.r + z$gabor.i * i: 2D array containing the gabor
#           functions located on the ridge
#
#########################################################################
{

   
   gabor.r <- numeric(nbnodes * sigsize)
   gabor.i <- numeric(nbnodes * sigsize)


   z <- .C("vgabor_time",
           as.double(frequency),
           as.double(scale),
           as.integer(location),
           gabor.r = as.double(gabor.r),
           gabor.i = as.double(gabor.i),
           as.integer(sigsize),
           as.integer(nbnodes),
           PACKAGE="Rwave")

   z$gabor.r + 1i*z$gabor.i
}











