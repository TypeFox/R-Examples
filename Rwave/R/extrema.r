#########################################################################
#    $log: extrema.S,v $
#########################################################################
#    (c) Copyright 1997
#               by
#   Author: Rene Carmona, Bruno Torresani, Wen L. Hwang, Andrea Wang
#              Princeton University 
#              All right reserved
########################################################################




ext <- function(wt, scale=FALSE, plot=TRUE)
#****************************************************************
# ext
# ---
#   compute the extrema of the dyadic wavelet transform.
#
# input
# -----
#   wt: wavelet transform
#   scale: a flag indicating if the extrema at each 
#	   resolution will be plotted at the same sacle.
#
# output
# ------
#   original: original signal
#   extrema: extrema representation
#   Sf: coarse resolution of signal
#   maxresoln: number of decomposition
#   np: size of signal
#****************************************************************
{
  s <- wt$original
  maxresoln <- wt$maxresoln
  np <- wt$np

  extrema <- matrix(0, nrow=maxresoln, ncol=np)
  extrema <- t(extrema)
  dim(extrema) <- c(length(extrema), 1)

  t(wt$Wf)   # Note: transposed wt is not assigned to wt
  dim(wt$Wf) <- c(length(wt$Wf), 1)

  z <- .C("modulus_maxima", 
          a=as.double(extrema),
          as.double(wt$Wf),
          as.integer(maxresoln), 
          as.integer(np),
           PACKAGE="Rwave")

  extrema <- t(z$a)
  dim(extrema) <- c(np, maxresoln)
  
  cat("number of extrema =", sum(extrema!=0), "\n")

  if(plot)  plotResult(extrema, s, maxresoln, scale)

    list(original=s,extrema=extrema,Sf=wt$Sf,maxresoln=maxresoln,np=np)
}












