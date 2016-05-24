#########################################################################
#    (c) Copyright 1997
#               by
#   Author: Rene Carmona, Bruno Torresani, Wen L. Hwang, Andrea Wang
#              Princeton University 
#              All right reserved
########################################################################

mrecons <- function(extrema, filtername="Gaussian1", readflag=FALSE)
#****************************************************************
#  mrecons
#  --------
#    Reconstruct from dyadic extrema. The reconstructed signal
#    preserves locations and values at extrema. This is Carmona's
#    method of extrema reconstruction.	
# 
#  input
#  -----
#    extrema : the extrema representation 
#    filtername : filters (Gaussain1 stands filters corresponding to
#	Mallat and Zhong's wavelet. And, the Haar stands for the
#       filters for Haar basis.	
#    readflag: if set to TRUE, read kernel from a precomputed file.
#
#  output
#  ------
#    f: signal reconstructed from the wavelet transform extrema
#    g: f plus mean of f
#    h: f plus coarse scale of f
#****************************************************************
{
  pure <- extrema$original
  x <- adjust.length(pure)
  pure.data <- x$signal

  data <- extrema$original
  maxresoln <- extrema$maxresoln
  np <- extrema$np

  if(2^maxresoln > np)
    stop("The support of the coarsest scale filter is bigger than signal size")
  if((readflag && (maxresoln > 9)) || (np > 4096))
    stop("Can not use readflag")

  dim(extrema$extrema) <- c(length(extrema$extrema), 1)

  f <- 1:np

  z <- .C("extrema_reconst",
          as.character(filtername),
          a = as.double(f),	
          as.double(extrema$extrema),
          as.integer(maxresoln),
          as.integer(np),
          as.integer(readflag),
           PACKAGE="Rwave")

  f <- z$a 
  g <- f + mean(pure.data)
  h <- f + extrema$Sf

  lim <- range(f, g, h, pure.data)
  par(mfrow=c(3,1))

  plot.ts(f, ylim=lim, xlab="", ylab="", main="Reconstruction")
  plot.ts(g, ylim=lim, xlab="", ylab="", main="Reconstruction + Mean")
  plot.ts(h, ylim=lim, xlab="", ylab="", main="Reconstruction + Sf")

  par(mfrow=c(1,1))

  list(f=f, g=g, h=h)
}

dwinverse <- function(wt, filtername="Gaussian1")
#****************************************************************
#  dwinverse
#  ---------
#    Invert the dyadic wavelet transform.
# 
#  input
#  -----
#    wt: the wavelet transform 
#    filtername : filters used. ("Gaussian1" stands for the filters
#       corresponds to those of Mallat and Zhong's wavlet. And
#       "Haar" stands for the filters of Haar basis.	
#
#  output
#  ------
#    f: the reconstructed signal
#****************************************************************
{
  Sf <- wt$Sf
  Sf <- c(Sf)	
  dim(Sf) <- c(length(Sf), 1)
  Wf <- wt$Wf
  Wf <- c(Wf)	
  dim(Wf) <- c(length(Wf), 1)
  np <- wt$np
  maxresoln <- wt$maxresoln
  fback <- 1:np
  dim(fback) <- c(length(fback), 1)

  z <- .C("Sinverse_wavelet_transform",
          a = as.double(fback),
          as.double(Sf),
          as.double(Wf),
          as.integer(maxresoln),
          as.integer(np), 
          as.character(filtername),
           PACKAGE="Rwave")

  f <- z$a 
  f
}




