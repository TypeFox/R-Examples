#########################################################################
#      $Log: WV.S,v $
#########################################################################
#
#               (c) Copyright  1998
#                          by                                   
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  Princeton University 
#                  All right reserved                           
#########################################################################





WV <- function(input, nvoice, freqstep = (1/nvoice), plot = TRUE)
#########################################################################
#       WV:   
#       ---
#        Wigner-Ville function
# 	  compute the Wigner-Ville transform, without any smoothing
#
#       input:
#       ------
# 	 input: input signal (possibly complex-valued)
#	 nvoice: number of frequency bands
#        freqstep: sampling rate for the frequency axis
#	 plot: if set to TRUE, displays the modulus of cwt on the graphic
#		device.
#
#       output:
#       -------
#        output: (complex) Wigner-Ville transform
#
#########################################################################
{
  oldinput <- input
  isize <- length(oldinput)
  
  tmp <- adjust.length(oldinput)
  input <- tmp$signal
  newsize <- length(input)
  
  ## pp <- nvoice
  pp <- newsize
  Routput <- matrix(0, newsize, pp)
  Ioutput <- matrix(0, newsize, pp)
  output <- matrix(0, newsize, pp)
  dim(Routput) <- c(pp * newsize, 1)
  dim(Ioutput) <- c(pp * newsize, 1)
  dim(input) <- c(newsize, 1)
  
  z <- .C("WV",
          as.double(input),
          Rtmp = as.double(Routput),
          Itmp = as.double(Ioutput),
          as.integer(nvoice),
          as.double(freqstep),
          as.integer(newsize),
          PACKAGE="Rwave")
  
  Routput <- z$Rtmp
  Ioutput <- z$Itmp
  dim(Routput) <- c(newsize, pp)
  dim(Ioutput) <- c(newsize, pp)
  
  output <- Routput[1:isize,] + 1i*Ioutput[1:isize,]
  if(plot) {
    image(Mod(output[,1:(isize/2)]),
          xlab="Time", ylab="Frequency")
    title("Wigner-Ville Transform Modulus")
  }
  output
}

