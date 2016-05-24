#########################################################################
#    $log: mw.S,v $
#########################################################################
#
#    (c) Copyright 1997
#               by
#   Author: Rene Carmona, Bruno Torresani, Wen L. Hwang, Andrea Wang
#              Princeton University 
#              All right reserved
########################################################################




check.maxresoln <- function( maxresoln, np )
#*********************************************************************#
# check.maxresoln
# ---------------
# stop when the size of 2^maxresoln is no less than signal size
#
# input 
# -----
# maxresoln: number of decomposition
# np: signal size
#
# output
# ------
#*********************************************************************#
{
  if ( 2^(maxresoln+1) > np )
    stop("maxresoln is too large for the given signal")
}


mw <- function(inputdata,maxresoln,filtername="Gaussian1",scale=FALSE,
	       plot=TRUE)
#*********************************************************************#
# mw
# --
#   mw computes the wavelet decomposition by Mallat and Zhong's wavelet.
#
# input
# -----
#   inputdata: either a text file or an S object of input data.
#   maxresoln: number of decomposition
#   filename: name of filter (Gaussian1 stands for Mallat and Zhong's filter ;
#             Haar stands for the Haar basis).	
#   scale: when set, the wavelet transform at each scale will be plotted 
#          with the same scale.
#   plot : indicate if the wavelet transform at each scale will be plotted.
#          		
# output
# ------
#   original: original signal
#   Wf: wavelet transform of signal
#   Sf: signal at a lower resolution
#   maxresoln: number of decomposition
#   np: size of signal
#*********************************************************************#
{
  x <- adjust.length(inputdata)
  s <- x$signal
  np <- x$length

  Sf <- matrix(0, nrow=(maxresoln+1), ncol=np)
  Wf <- matrix(0, nrow=maxresoln, ncol=np)  

  ## Convert a matrix maxresoln by np matrix into a vector
  Sf <- t(Sf)
  dim(Sf) <- c(length(Sf),1)
  Wf <- t(Wf)
  dim(Wf) <- c(length(Wf),1)

  y <- .C("Sf_compute", 
	Sf=as.double(Sf), 
	as.double(s), 
	as.integer(maxresoln), 
	as.integer(np),
	as.character(filtername),
           PACKAGE="Rwave")

  z <- .C("Wf_compute",
	Wf=as.double(Wf),
	as.double(y$Sf),
	as.integer(maxresoln),
	as.integer(np),
	as.character(filtername),
           PACKAGE="Rwave")

  ## Convert the vectors into the original matrix
  Sf <- t(y$Sf)

  Sf <- Sf[(maxresoln*np+1):((maxresoln+1)*np)]

  Wf <- t(z$Wf)
  dim(Wf) <- c(np, maxresoln)

  if(plot)
    plotwt(s, Wf, Sf, maxresoln, scale)

  list( original=s, Wf=Wf, Sf=Sf, maxresoln=maxresoln, np=np )
}







