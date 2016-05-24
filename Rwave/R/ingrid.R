#########################################################################
#
#               (c) Copyright  1997
#                          by                                   
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  Princeton University
#                  All right reserved                           
#########################################################################

dw <- function(inputdata, maxresoln, scale=FALSE, NW=6, plot=TRUE)
#*********************************************************************#
# dw computes the Daubechies wavelet transform from resolution 1 to
# the given maximum resolution.  If phi.flag is TRUE, it returns a 
# list consisting psi and phi wavelet transforms as two matrices.
#*********************************************************************#
{
  if(! ((NW > 1) && (NW <11)))
    stop("NW has to be between 2 and 10")
  
  x <- adjust.length(inputdata)
  s <- x$signal
  np <- x$length
  check.maxresoln(maxresoln, np)
  phi <- matrix(0, nrow=maxresoln+1, ncol=np)
  psi <- matrix(0, nrow=maxresoln, ncol=np)

  ## Convert matrices into vectors
  phi <- t(phi)
  dim(phi) <- c(length(phi), 1)
  psi <- t(psi)
  dim(psi) <- c(length(psi), 1)
  
  z <- .C("daubechies_wt", 
          phi=as.double(phi),
          psi=as.double(psi),
          as.double(s),
          as.integer(NW),
          as.integer(maxresoln),
          as.integer(np),
           PACKAGE="Rwave")
  
  ## Convert vectors into original matrices
  phi <- t(z$phi)
  dim(phi) <- c(np, maxresoln+1)
  psi <- t(z$psi)
  dim(psi) <- c(np, maxresoln)
  
  if(plot)
    plotwt(s, psi, phi, maxresoln, scale)

  phi <- matrix(phi[,(maxresoln+1)], ncol=1, nrow=np)
  list(original=s, Wf=psi, Sf=phi, maxresoln=maxresoln, np=np, NW=NW)
}

ddw <- function(inputdata, maxresoln, scale=FALSE, NW=6)
#*********************************************************************#
# ddw computes the discrete Daubechies wavelet transform from 
# resolution 1 to the given maximum resolution.  
#*********************************************************************#
{
  if(! ((NW > 1) && (NW < 11)))
    stop("NW has to be between 2 and 10")

  x <- adjust.length(inputdata)
  s <- x$signal
  np <- x$length	
  check.maxresoln(maxresoln, np)
  phi <- matrix(0, nrow=maxresoln+1, ncol=np)
  psi <- matrix(0, nrow=maxresoln, ncol=np)

  ## Convert matrices into vectors
  phi <- t(phi)
  dim(phi) <- c(length(phi), 1)
  psi <- t(psi)
  dim(psi) <- c(length(psi), 1)
  
  z <- .C("compute_ddwave", 
          phi=as.double(phi),
          psi=as.double(psi),
          as.double(s),
          as.integer(maxresoln),
          as.integer(np),
          as.integer(NW),
           PACKAGE="Rwave")
  
  ## Convert vectors into original matrices
  phi <- t(z$phi)
  dim(phi) <- c(np, maxresoln+1)
  psi <- t(z$psi)
  dim(psi) <- c(np, maxresoln)
  plotwt(s, psi, phi, maxresoln, scale)
  phi <- matrix(phi[,(maxresoln+1)], ncol=1, nrow=np)
  list(original=s, Wf=psi, Sf=phi, maxresoln=maxresoln, np=np, NW=NW)
}
