#########################################################################
#      $Log: Cwt_Morlet.S,v $
#########################################################################
#
#               (c) Copyright  1997                             
#                          by                                   
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  Princeton University
#                  All right reserved                           
#########################################################################

cwt <- function(input, noctave, nvoice = 1, w0 = 2*pi, twoD = TRUE,
                plot = TRUE)
#########################################################################
#      cwt:
#      ---
#       continuous wavelet transform function
# 	compute the continuous wavelet transform with (complex-valued)
#			Morlet wavelet
#
#       Input:
#       ------
# 	 input: input signal (possibly complex-valued)
#	 noctave: number of powers of 2 for the scale variable
#	 nvoice: number of scales between 2 consecutive powers of 2
#        w0: central frequency of Morlet wavelet
#	 twoD: if set to TRUE, organizes the output as a 2D array 
#			(signal_size X nb_scales)
#		      if not: 3D array (signal_size X noctave X nvoice)
#	 plot: if set to TRUE, displays the modulus of cwt on the graphic
#		device.
#
#       Output:
#       -------
#        tmp: continuous (complex) wavelet transform
#
#########################################################################
{
  oldinput <- input
  isize <- length(oldinput)
  
  tmp <- adjust.length(oldinput)
  input <- tmp$signal
  newsize <- length(input)
  
  pp <- noctave * nvoice
  Routput <- matrix(0,newsize,pp)
  Ioutput <- matrix(0,newsize,pp)
  output <- matrix(0,newsize,pp)
  dim(Routput) <- c(pp * newsize,1)
  dim(Ioutput) <- c(pp * newsize,1)
  dim(input) <- c(newsize,1)
  
  z <- .C("Scwt_morlet",
          as.double(Re(input)),
          as.double(Im(input)),
          Rtmp = as.double(Routput),
          Itmp = as.double(Ioutput),
          as.integer(noctave),
          as.integer(nvoice),
          as.integer(newsize),
          as.double(w0),
          PACKAGE="Rwave")
  
  Routput <- z$Rtmp
  Ioutput <- z$Itmp
  dim(Routput) <- c(newsize,pp)
  dim(Ioutput) <- c(newsize,pp)
  if(twoD) {
    output <- Routput[1:isize,] + 1i*Ioutput[1:isize,]
    if(plot) {
      image(Mod(output), xlab="Time", ylab="log(scale)",
            main="Wavelet Transform Modulus")
    }
    output
  } 
  else {
    Rtmp <- array(0,c(isize,noctave,nvoice))
    Itmp <- array(0,c(isize,noctave,nvoice))
    for(i in 1:noctave)
      for(j in 1:nvoice) {
        Rtmp[,i,j] <- Routput[1:isize,(i-1)*nvoice+j]
        Itmp[,i,j] <- Ioutput[1:isize,(i-1)*nvoice+j]
      }
    Rtmp + 1i*Itmp
  }
}

cwtpolar <- function(cwt, threshold=0.0)
#########################################################################
#       cwtpolar:   
#       --------
# 	continuous wavelet transform conversion:
#        converts one of the possible outputs of cwt to modulus and phase
#
#       input:
#       ------
# 	 cwt: 3D array containing a continuous wavelet transform (output
#		of cwt, with twoD=FALSE)
#	 threshold: the phase is forced to -pi if the modulus is 
#		less than threshold.
#
#       output:
#       -------
#        output1: modulus
#        output2: phase      
#
#########################################################################
{
  tmp1 <- cwt
  sigsize <- dim(tmp1)[1] # sig size
  noctave <- dim(tmp1)[2]
  nvoice <- dim(tmp1)[3]
  
  output1 <- array(0,c(sigsize,noctave,nvoice))
  output2 <- array(0,c(sigsize,noctave,nvoice))
  for(i in 1:noctave)
    for(j in 1:nvoice) {
      output1[,i,j] <- sqrt(Re(tmp1[,i,j])^2 + Im(tmp1[,i,j])^2)
      output2[,i,j] <- atan2(Im(tmp1[,i,j]), Re(tmp1[,i,j]))
    }
  
  ma <- max(output1)
  rel <- threshold * ma
  dim(output1) <- c(sigsize * noctave * nvoice,1)
  dim(output2) <- c(sigsize * noctave * nvoice,1)
  output2[abs(output1) < rel] <- -3.14159
  dim(output1) <- c(sigsize,noctave,nvoice)
  dim(output2) <- c(sigsize,noctave,nvoice)
  
  list(modulus=output1,argument=output2)
}
 
cwtimage <- function(input)
#########################################################################
#       cwtimage:   
#       ---------
#        continuous wavelet transform display
# 	 converts the output (modulus or argument) of cwtpolar to a 
#		2D array and displays on the graphic device
#
#       input:
#       ------
# 	 input: 3D array containing a continuous wavelet transform (output
#		of cwtpolar
#       output:
#       -------
#        output: 2D array 
#
#########################################################################
{
  sigsize <- dim(input)[1]
  noctave <- dim(input)[2]
  nvoice <- dim(input)[3]
  
  output <- matrix(0,sigsize, noctave * nvoice)
  for(i in 1:noctave) {
    k <- (i-1) * nvoice
    for(j in 1:nvoice) {
      output[,k+j] <- input[,i,j]
    }
  }
  image(output)
  output
}
 
vwt <- function(input, scale, w0 = 2*pi)
#########################################################################
#      vwt:   voice Morlet wavelet transform   
#      ----
#       continuous wavelet transform on one scale
# 	compute the continuous wavelet transform with (complex-valued)
#			Morlet wavelet
#
#       input:
#       ------
# 	 input: input signal (possibly complex-valued)
#        scale: value of the scale at which the transform is computed
#        w0: central frequency of Morlet wavelet
#
#       output:
#       -------
#        Routput + i Ioutput: voice wavelet transform (complex 1D array)
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
  
  z <- .C("Svwt_morlet",
          as.double(Re(input)),
          as.double(Im(input)),
          Rtmp = as.double(Routput),
          Itmp = as.double(Ioutput),
          as.double(scale),
          as.integer(newsize),
          as.double(w0),
          PACKAGE="Rwave")
  
  Routput <- z$Rtmp
  Ioutput <- z$Itmp
  Routput[1:isize] + 1i*Ioutput[1:isize]
}

morlet <- function(sigsize, location, scale, w0 = 2*pi)
#########################################################################
#       morlet: Morlet's wavelet  
#       -------
#        Generates a Morlet wavelet for given location and scale
#
#       input:
#       ------
#        sigsize: signal size (dimension of the array)
#        location: location of the wavelet
#        scale: value of the scale at which the transform is computed
#        w0: central frequency of Morlet wavelet
#
#       output:
#       -------
#        z$wavelet.r + z$wavelet.i * i: wavelet (complex 1D array
#          of size sigsize)
#
#########################################################################
{
  wavelet.r <- numeric(sigsize)
  wavelet.i <- numeric(sigsize)
  
  z <- .C("morlet_time",
          as.double(w0),
          as.double(scale),
          as.integer(location),
          wavelet.r = as.double(wavelet.r),
          wavelet.i = as.double(wavelet.i),
          as.integer(sigsize),
          PACKAGE="Rwave")
  
  z$wavelet.r + 1i*z$wavelet.i
}

vecmorlet <- function(sigsize, nbnodes, bridge, aridge, w0 = 2*pi)
#########################################################################
#       vecmorlet:   vector of Morlet's wavelets
#       ---------
#        Generates morlet wavelets located on the ridge
#
#       input:
#       ------
#        sigsize: signal size (dimension of the array)
#        nbnodes: number of ridge samples
#        bridge: b coordinates of the ridge samples
#        aridge: acoordinates of the ridge samples
#        w0: central frequency of Morlet wavelet
#
#       output:
#       -------
#        z$morlet.r + z$morlet.i * i: 2D array containing the morlet
#           wavelets located on the ridge
#
#########################################################################
{
  morlet.r <- numeric(nbnodes * sigsize)
  morlet.i <- numeric(nbnodes * sigsize)

  z <- .C("vmorlet_time",
          as.double(w0),
          as.double(aridge),
          as.integer(bridge),
          morlet.r = as.double(morlet.r),
          morlet.i = as.double(morlet.i),
          as.integer(sigsize),
          as.integer(nbnodes),
          PACKAGE="Rwave")
  
  z$morlet.r + 1i*z$morlet.i
}
