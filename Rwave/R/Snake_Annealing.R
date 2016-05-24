#########################################################################
#      $Log: Snake_Annealing.S,v $
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



snake <- function(tfrep, guessA, guessB, snakesize = length(guessB),
                  tfspec = numeric(dim(modulus)[2]), subrate = 1,
                  temprate = 3, muA = 1, muB = muA, lambdaB = 2*muB,
                  lambdaA = 2*muA, iteration = 1000000, seed = -7,
                  costsub = 1, stagnant = 20000, plot = TRUE)
#########################################################################
#     snake:   
#     ------
#      ridge extraction using snakes and simulated annealing
#
#       input:
#       ------
# 	 tfrep: the wavelet or Gabor transform.
#	 guessA: initial guess for ridge function (frequency coordinate).
#	 guessB: initial guess for ridge function (time coordinate).
#        tfspec: estimate for the contribution of the noise to the
#               transform modulus (1D array)
#        subrate: subsampling rate for the transform modulus.
#	 temprate: constant (numerator) in the temperature schedule.
#	 muA, muB: coeff of phi' and rho' in cost function
#	 lambdaA,lambdaB: same thing for 2nd derivatives
#        iteration: maximal number of iterations
#        seed: initialization for random numbers generator
#        costsub: subsampling of the cost function in cost
#               costsub=1 means that the whole cost function
#               is returned
#        stagnant: maximum number of steps without move (for the
#                  stopping criterion)
#        plot :	 when set (by default), certain results will be
#	         displayed	
#
#       output:
#       ------
#        A: 1D array containing the frequency coordinate of the ridge
#        B: 1D array containing the time coordinate of the ridge
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
  phi <- as.integer(guessA)
  rho <- as.integer(guessB)
  count <- 0	
  dim(phi) <- c(length(phi),1)	
  dim(rho) <- c(length(rho),1)	
  if(plot== TRUE)  image(modulus)
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


   z <- .C("Ssnake_annealing",
        cost = as.double(cost),
        as.double(smodulus),
        phi = as.double(phi),
        rho = as.double(rho),
        as.double(lambdaA),
	as.double(muA),
        as.double(lambdaB),
	as.double(muB),
        as.double(temprate),
	as.integer(sigsize),
	as.integer(snakesize),
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
  cat("Number of moves:",count,"(",stagnant," still steps)\n")
  if(plot==TRUE) lines(z$rho+1, z$phi+1)
  list(A = z$phi+1, B = z$rho+1, cost = z$cost[1:count])

}
  


   
snakoid <- function(modulus, guessA, guessB, snakesize = length(guessB),
                    tfspec = numeric(dim(modulus)[2]), subrate = 1,
                    temprate = 3, muA = 1, muB = muA, lambdaB = 2*muB,
                    lambdaA = 2*muA, iteration = 1000000, seed = -7,
                    costsub = 1, stagnant = 20000, plot = TRUE)
#########################################################################
#     snakoid:   
#     ----------
#      ridge extraction using snakes and simulated annealing
#
#       input:
#       ------
# 	 modulus: modulus of the wavelet or Gabor transform.
#	 guessA: initial guess for ridge function (frequency coordinate).
#	 guessB: initial guess for ridge function (time coordinate).
#        tfspec: estimate for the contribution of the noise to the
#               transform modulus (1D array)
#        subrate: subsampling rate for the transform modulus.
#	 temprate: constant (numerator) in the temperature schedule.
#	 muA, muB: coeff of phi' and rho' in cost function
#	 lambdaA,lambdaB: same thing for 2nd derivatives
#        iteration: maximal number of iterations
#        seed: initialization for random numbers generator
#        costsub: subsampling of the cost function in cost
#               costsub=1 means that the whole cost function
#               is returned
#        stagnant: maximum number of steps without move (for the
#                  stopping criterion)
#        plot:	 when set(default), some results will be displayed	
#
#       output:
#       ------
#        A: 1D array containing the frequency coordinate of the ridge
#        B: 1D array containing the time coordinate of the ridge
#        cost: 1D array containing the cost function
#
#########################################################################
{

  sigsize <- dim(modulus)[1]
  nscale <- dim(modulus)[2]
  costsize <- as.integer(iteration/costsub) + 1	
  cost <- 1:costsize
  tfspectrum <- tfspec 

  cost[] <- 0
  phi <- guessA
  rho <- guessB
  count <- 0	
  dim(phi) <- c(length(phi),1)	
  dim(rho) <- c(length(rho),1)
  if(plot== TRUE)  image(modulus)	
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
        as.double(subrate),
           PACKAGE="Rwave")

  smodulus <- z$smodulus
  smodsize <- z$smodsize
  smodulus <- smodulus * smodulus
  dim(smodulus) <- c(smodsize, nscale)
  for (k in 1:nscale)
    smodulus[,k] <- smodulus[,k] - tfspectrum[k]
  dim(smodulus) <- c(smodsize * nscale, 1)


   z <- .C("Ssnakenoid_annealing",
        cost = as.double(cost),
        as.double(smodulus),
        phi = as.double(phi),
        rho = as.double(rho),
        as.double(lambdaA),
	as.double(muA),
        as.double(lambdaB),
	as.double(muB),
        as.double(temprate),
	as.integer(sigsize),
	as.integer(snakesize),
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
  cat("Number of moves:",count,"(",stagnant," still steps)\n")
  if(plot== TRUE)  lines(z$rho+1, z$phi+1) 
  list(A = z$phi+1, B = z$rho+1, cost = z$cost[1:count])

}
  


snakeview <- function(modulus,snake)
#########################################################################
#      snakeview:
#      ----------
#       Restrict time-frequency transform to a snake
#
#       input:
#       ------
#	 modulus: modulus of the wavelet (or Gabor) transform
#        snake: B and A components of a snake
#
#       output:
#       -------
#        img: 2D array containing the restriction of the transform
#             modulus to the snake
#
#########################################################################
{
   A <- snake$A
   B <- snake$B
   snakesize <- length(A)
   sigsize <- dim(modulus)[1]
   nscale <- dim(modulus)[2]

   img <- matrix(0,sigsize,nscale)
   for(i in 1:snakesize)
     img[B[i],A[i]] <- modulus[B[i],A[i]]

   img
}

     









