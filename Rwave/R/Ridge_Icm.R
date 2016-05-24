#########################################################################
#      $Log: Ridge_Icm.S,v $
#########################################################################
#
#               (c) Copyright  1997                             
#                          by                                   
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  Princeton University 
#                  All right reserved                           
#########################################################################




icm <- function(modulus, guess, tfspec = numeric(dim(modulus)[2]),
	subrate = 1, mu = 1, lambda = 2*mu, iteration = 100)
#########################################################################
#       icm:   
#       ----
#        ridge extraction using icm algorithm
#
#       Input:
#       ------
# 	 modulus: modulus of the wavelet or Gabor transform.
#	 guess: initial guess for ridge function (1D array).
#        tfspec: estimate for the contribution of the noise to the
#               transform square modulus (1D array).
#        subrate: subsampling rate for the transform square modulus.
#        lambda: coefficient in front of phi'' in the cost function
#        mu: coefficient in front of phi' in the cost function
#        iteration: maximal number of iterations
#
#       Output:
#       -------
#        ridge: 1D array (of length sigsize) containing the ridge
#        cost: 1D array containing the cost function
#
#########################################################################
{

  sigsize <- dim(modulus)[1]
  nscale <- dim(modulus)[2]
  costsize <- iteration	
  cost <- 1:costsize
  tfspectrum <- tfspec

  cost[] <- 0
  phi <- as.integer(guess)
  count <- 0	
  dim(phi) <- c(sigsize,1)	
  dim(modulus) <- c(sigsize * nscale, 1)

  smodsize <- as.integer(sigsize/subrate)

  if((sigsize/subrate - as.integer(sigsize/subrate)) > 0)
    smodsize <- as.integer(sigsize/subrate) + 1

  smodulus <- matrix(0,smodsize,nscale)
  dim(smodulus) <- c(smodsize * nscale, 1)

  if (subrate != 1){
  z <- .C("Smodulus_smoothing",
          as.double(modulus),
          smodulus = as.double(smodulus),
          as.integer(sigsize),
          smodsize = as.integer(smodsize),            
          as.integer(nscale),
          as.integer(subrate),
           PACKAGE="Rwave")
  
  smodulus <- z$smodulus
  smodsize <- z$smodsize
  }
  else smodulus <- modulus
  smodulus <- smodulus * smodulus
  dim(smodulus) <- c(smodsize, nscale)
  for (k in 1:nscale)
    smodulus[,k] <- smodulus[,k] - tfspectrum[k]
  dim(smodulus) <- c(smodsize * nscale, 1)

  z <- .C("Sridge_icm",
          cost = as.double(cost),
          as.double(smodulus),
          phi = as.double(phi),
          as.double(lambda),
          as.double(mu),
          as.integer(sigsize),
          as.integer(nscale),
          as.integer(iteration),
          nb = as.integer(count),
          as.integer(subrate),
          as.integer(smodsize),
           PACKAGE="Rwave")

  count <- z$nb
  cat("Number of iterations:",count,"\n")
  list(ridge = z$phi+1, cost = z$cost[1:count])
}
  




















