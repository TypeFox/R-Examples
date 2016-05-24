#########################################################################
#       $Log: TF_Maxima.S,v $
#
#               (c) Copyright  1997
#                          by                                   
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  Princeton University 
#                  All right reserved                           
#########################################################################





tfgmax <- function(input, plot = TRUE)
#########################################################################
#  tfgmax:  
#  -------
# 	 Continuous time-frequency transform global maxima:
#          compute the continuous wavelet transform global maxima (for
#    	   fixed position)
#
#      input:
#      ------
# 	input: continuous time-frequency transform (2D array)
#	plot: if set to TRUE, displays the maxima of cwt on the graphic
#          device.
#
#      output:
#      -------
#       output: values of the maxima (1D array)
#       pos: positions of the maxima (1D array)
#       
#########################################################################
{

  sigsize <- dim(input)[1]
  pp <- dim(input)[2]
  input1 <- input
	
  output <- matrix(0,sigsize,pp)
  dim(input1) <- c(sigsize * pp,1)
  dim(output) <- c(sigsize * pp,1)
  posvector <- 1:sigsize
  posvector[] <- 0

  z <- .C("Scwt_gmax",
           as.double(input1),
           output = as.double(output),
           as.integer(sigsize),
           as.integer(pp),
           pos = as.integer(posvector),
           PACKAGE="Rwave")

  output <- z$output
  pos <- z$pos

  dim(output) <- c(sigsize, pp)
  if(plot)image(output)
  list(output=output, pos= pos)
}



tflmax <- function(input, plot = TRUE)
#########################################################################
#  tflmax:
#  -------
# 	continuous time-frequency transform local maxima:
#        compute the time-frequency transform local maxima (for
#    	 fixed position)
#
#      input:
#      ------
# 	input: continuous time-frequency transform (2D array)
#	plot: if set to TRUE, displays the maxima of cwt on the graphic
#          device.
#
#      output:
#      -------
#       output: values of the maxima (2D array)
#
#########################################################################
{

  sigsize <- dim(input)[1]
  pp <- dim(input)[2]
  input1 <- input

  output <- matrix(0,sigsize,pp)
  dim(input1) <- c(sigsize * pp,1)
  dim(output) <- c(sigsize * pp,1)

  z <- .C("Scwt_mridge",
           as.double(input1),
           output = as.double(output),
           as.integer(sigsize),
           as.integer(pp),
           PACKAGE="Rwave")

  output <- z$output
  dim(output) <- c(sigsize, pp)
  if(plot)image(output)
  output
}



cleanph <- function(tfrep, thresh = .01, plot = TRUE)
#########################################################################
#  cleanph:
#  --------
# 	sets to zero the phase of time-frequency transform when
#        modulus is below a certain value.
#
#      input:
#      ------
# 	tfrep: continuous time-frequency transform (2D array)
#       thresh: (relative) threshold.
#	plot: if set to TRUE, displays the maxima of cwt on the graphic
#          device.
#
#      output:
#      -------
#       output: thresholded phase (2D array)
#
#########################################################################
{
  thrmod1 <- Mod(tfrep)
  thrmod2 <- Mod(tfrep)
  limit <- range(thrmod1)[2] * thresh
  thrmod1 <- (Mod(tfrep) > limit)
  thrmod2 <- (Mod(tfrep) <= limit)

  output <- thrmod1 * Arg(tfrep) - pi * thrmod2

  if(plot) image(output)
  output
}

























