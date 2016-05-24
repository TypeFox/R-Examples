#########################################################################
#       $Log: Noise.S,v $
# Revision 1.2  1995/04/05  18:56:55  bruno
# *** empty log message ***
#
# Revision 1.1  1995/04/02  01:04:16  bruno
# Initial revision
#
#               (c) Copyright  1997
#                          by                                   
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  Princeton University
#                  All right reserved                           
#########################################################################




tfpct <- function(input,percent=.8,plot=TRUE)
#########################################################################
#       tfpct:   
#       ------
# 	 compute a percentile of time-frequency representation frequency 
#		by frequency
#
#       input:
#       ------
# 	 input: modulus of the continuous wavelet transform (2D array)
#	 percent: value of the percentile to be computed
#        plot: if set to TRUE, displays the values of the energy as a 
#                  function of the scale
#
#       output:
#       -------
#        output: 1D array of size nbscales containing the noise estimate
#
#########################################################################
{

  nscale <- dim(input)[2]

  output <- numeric(nscale)
  for(i in 1:nscale) output[i] <- quantile(input[,i],percent)

  if(plot)plot.ts(output)
  output
}


tfmean <- function(input,plot=TRUE)
#########################################################################
#       tfmean:   
#       -------
# 	 compute the mean of time-frequency representation frequency 
#		by frequency
#
#       input:
#       ------
# 	 input: modulus of the continuous wavelet transform (2D array)
#        plot: if set to TRUE, displays the values of the energy as a 
#                  function of the scale
#
#       output:
#       -------
#        output: 1D array of size nbscales containing the noise estimate
#
#########################################################################
{

  nscale <- dim(input)[2]

  output <- numeric(nscale)
  for(i in 1:nscale) output[i] <- mean(input[,i])

  if(plot) plot.ts(output)
  output
}


tfvar <- function(input,plot=TRUE)
#########################################################################
#       tfvar:   
#       ------
# 	 compute the variance of time-frequency representation frequency 
#		by frequency
#
#       input:
#       ------
# 	 input: modulus of the continuous wavelet transform (2D array)
#        plot: if set to TRUE, displays the values of the energy as a 
#                  function of the scale
#
#       output:
#       -------
#        output: 1D array of size nbscales containing the noise estimate
#
#########################################################################
{

  nscale <- dim(input)[2]

  output <- numeric(nscale)
  for(i in 1:nscale) output[i] <- var(input[,i])

  if(plot) plot.ts(output)
  output
}




















