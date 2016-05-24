#########################################################################
#      $Log: Ridge_Annealing.S,v $
#########################################################################
#
#               (c) Copyright  1997
#                          by                                   
#      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   
#                  University of California, Irvine             
#                  All right reserved                           
#########################################################################





corona <- function(tfrep, guess, tfspec = numeric(dim(tfrep)[2]),
                   subrate = 1, temprate = 3, mu = 1, lambda = 2*mu,
                   iteration = 1000000, seed = -7, stagnant = 20000,
                   costsub = 1, plot = TRUE)
#########################################################################
#       corona:
#       -------
#        ridge extraction using simulated annealing
# 	  compute the continuous wavelet transform ridge from the cwt
#	   modulus, using simulated annealing
#
#       Input:
#       ------
# 	 tfrep: wavelet or Gabor transform.
#	 guess: initial guess for ridge function (1D array).
#        tfspec: estimate for the contribution of the noise to the
#               transform modulus (1D array)
#        subrate: subsampling rate for the transform modulus.
#	 temprate: constant (numerator) in the temperature schedule.
#        lambda: coefficient in front of phi' in the cost function
#        mu: coefficient in front of phi'' in the cost function
#        iteration: maximal number of iterations
#        seed: initialization for random numbers generator
#        stagnant: maximum number of steps without move (for the
#                  stopping criterion)
#        costsub: subsampling of the cost function in cost
#               costsub = 1 means that the whole cost function
#               is returned
#        plot:	 when set(default), some results will be shown on 
#	        the display	
#
#       Output:
#       -------
#        ridge: 1D array (of length sigsize) containing the ridge
#        cost: 1D array containing the cost function
#
#########################################################################
{

  sigsize <- dim(tfrep)[1]
  nscale <- dim(tfrep)[2]
  blocksize <- 1	
  costsize <- as.integer(iteration/blocksize) + 1	
  cost <- 1:costsize

  modulus <- Mod(tfrep)
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
  smodulus <- smodulus * smodulus
  dim(smodulus) <- c(smodsize, nscale)
  for (k in 1:nscale)
    smodulus[,k] <- smodulus[,k] - tfspectrum[k]

  dim(smodulus) <- c(smodsize * nscale, 1)
  cat("dimension of smodulus",dim(smodulus),"\n")

   z <- .C("Sridge_annealing",
           cost = as.double(cost),
           as.double(smodulus),
           phi = as.double(phi),
           as.double(lambda),
           as.double(mu),
           as.double(temprate),
           as.integer(sigsize),
           as.integer(nscale),
           as.integer(iteration),
           as.integer(stagnant),
           as.integer(seed),
           nb = as.integer(count),
           as.integer(subrate),
           as.integer(costsub),
           as.integer(smodsize),
           PACKAGE="Rwave")

  count <- z$nb
  cat("Number of iterations:",(count-1)*costsub,"(stagnant=",stagnant,")\n")
  if(plot == TRUE) {
   image(Mod(tfrep))
   lines(z$phi+1)
  }	
  list(ridge = z$phi+1, cost = z$cost[1:count])
}



  
coronoid <- function(tfrep, guess, tfspec = numeric(dim(tfrep)[2]),
	subrate = 1, temprate = 3, mu = 1, lambda = 2*mu,
	iteration = 1000000, seed = -7, stagnant = 20000,
	costsub = 1,plot=TRUE)
#########################################################################
#       coronoid:   
#       ----------
#        ridge extraction using modified simulated annealing
#	   Modification of Sridge_annealing, the cost function is
#	   replaced with another one, in which the smoothness penalty
#	   is proportional to the transform modulus.
#
#       Input:
#       ------
# 	 tfrep: wavelet or Gabor transform.
#	 guess: initial guess for ridge function (1D array).
#        tfspec: estimate for the contribution of the noise to the
#               transform modulus (1D array)
#        subrate: subsampling rate for the transform modulus.
#	 temprate: constant (numerator) in the temperature schedule.
#        lambda: coefficient in front of phi' in the cost function
#        mu: coefficient in front of phi'' in the cost function
#        iteration: maximal number of iterations
#        seed: initialization for random numbers generator
#        stagnant: maximum number of steps without move (for the
#                  stopping criterion)
#        costsub: subsampling of the cost function in cost
#               costsub=1 means that the whole cost function
#               is returned
#        plot:	 when set(default), some results will be shown on 
#	        the display	
#
#       Output:
#       -------
#        ridge: 1D array (of length sigsize) containing the ridge
#        cost: 1D array containing the cost function
#
#########################################################################
{

  sigsize <- dim(tfrep)[1]
  nscale <- dim(tfrep)[2]
  costsize <- as.integer(iteration/costsub) + 1	
  cost <- 1:costsize

  modulus <- Mod(tfrep)
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
  smodulus <- smodulus * smodulus
  dim(smodulus) <- c(smodsize, nscale)
  for (k in 1:nscale)
    smodulus[,k] <- smodulus[,k] - tfspectrum[k]
  dim(smodulus) <- c(smodsize * nscale, 1)


   z <- .C("Sridge_coronoid",
        cost = as.double(cost),
        as.double(smodulus),
        phi = as.double(phi),
        as.double(lambda),
	as.double(mu),
        as.double(temprate),
	as.integer(sigsize),
        as.integer(nscale),
        as.integer(iteration),
	as.integer(stagnant),
	as.integer(seed),
	nb = as.integer(count),
        as.integer(subrate),
	as.integer(costsub),
        as.integer(smodsize),
           PACKAGE="Rwave")

  count <- z$nb
  cat("Number of iterations:",(count-1)*costsub,"(stagnant=",stagnant,")\n")
  if(plot == TRUE) {	
    image(Mod(tfrep))
    lines(z$phi+1)
  }
  list(ridge = z$phi+1, cost = z$cost[1:count])
}
  























