# copyright (C) 2014-2016 A.Rebecq and M.Chevalier
#### All functions in this file are private methods used by the
#### "calibration" function

# Function Wrapper
calib <- function(Xs, d, total, q=NULL, method=NULL, bounds = NULL,
                  alpha = NULL,
                  maxIter=500, calibTolerance=1e-06) {

  if(!is.null(method)) {
  switch(method,
          linear={
            inverseDistance <- inverseDistanceLinear
            params <- NULL
            updateParameters <- identity
          },
          raking={
            inverseDistance <- inverseDistanceRaking
            params <- NULL
            updateParameters <- identity
          },
          logit={
            inverseDistance <- inverseDistanceLogit
            params <- bounds
            updateParameters <- identity
          },
#            truncated={
#              inverseDistance <- inverseDistanceTruncated
#              params <- bounds
#              updateParameters <- identity
#            },
#           curlingHat={
#             inverseDistance <- distanceCurlingHat
#             # TODO : check params in list are correctly entered
#             params <- c(0.5,1.5) # For tests only
#             updateParameters <- updateParametersCurlingHat
#           },
          {
            print('By default, raking method selected')
            params <- NULL
            inverseDistance <- inverseDistanceRaking
            updateParameters <- identity
          }
  )
  } else {
    print('By default, raking method selected')
    params <- NULL
    inverseDistance <- inverseDistanceRaking
    updateParameters <- identity
  }
  # TODO : additional checks ?

  return(calibAlgorithm(Xs, d, total, q, inverseDistance,
                        updateParameters, params, maxIter, calibTolerance))

}

calibAlgorithm <- function(Xs, d, total, q=NULL,
                            inverseDistance, updateParameters, params,
                            maxIter=500, calibTolerance=1e-06) {

  if(is.null(q)) {
    q <- rep(1,length(d))
  }

  g <- NULL
  toleranceGInv = .Machine$double.eps # Tolerance when we compute ginv

  ## Linear optimization algorithm
  lambda = as.matrix(rep(0, ncol(Xs)))
  wTemp = as.vector(d * inverseDistance(Xs %*% lambda * q, params)[[1]])
  wUpdate <- d

  cont <- TRUE
  l <- 1

  while (cont) {

    phi = t(Xs) %*% wTemp - total
    T1 = t(Xs * wUpdate)
    phiprim = T1 %*% Xs
   
    wTemp <- NA
    
    try({
      lambda = lambda - ginv(phiprim, tol = toleranceGInv) %*% phi
      wTemp = as.vector(d * inverseDistance(Xs %*% lambda * q, params)[[1]])
      }
    )
	
	wUpdate <- as.vector(d * inverseDistance(Xs %*% lambda * q, params)[[2]])

    if (any(is.na(wTemp)) | any(is.infinite(wTemp))) {
      warning("No convergence")
      return(NULL)
    }

    tHat = t(Xs) %*% wTemp
    if (max(abs(tHat - total)/total) < calibTolerance) {
      cont <- FALSE
    }

    if(l >= maxIter) {
      cont <- FALSE
      warning(paste("No convergence in ", maxIter, " iterations."))
  	  return(NULL)
  	  return(wTemp)
    }

    l <- l+1
    # update Parameters before going back into loop
    updateParameters(params)

  }

  ## Return solution if found
  if(l <= maxIter) {
  	g = wTemp/d
    return(g)
  } else {
    warning(paste("No convergence in ", maxIter, " iterations."))
    return(NULL)
  }

  return(g)
}


# TODO: add qk vectors
# TODO: Explain why params in docs (uses S3 and not S4 methods)

##### Inverse distance functions
##### Each of these functions returns a list of the closed
##### form of the inverse functions of the distance used AND
##### its derivative.
inverseDistanceLinear <- function(x, params=NULL) {
  return(list((1+x),0))
}

inverseDistanceRaking <- function(x, params=NULL) {
  return(list(exp(x),exp(x)))
}

# Params are bounds
inverseDistanceLogit <- function(x, bounds) {

  if(length(bounds) != 2) {
    stop("Must enter LO and UP bounds in a vector.
          Example : bounds=c(0.5,1.5) for LO = 0.5 and UP=1.5")
  }

  L = bounds[1]
  U = bounds[2]
  A = (U - L) / ((1-L)*(U-1))
  distance = ( L*(U-1) + U*(1-L)*exp(A*x) ) / ( (U-1) + (1-L)*exp(A*x))

  fprim <- ( ( (U-L)**2 )*exp(A*x) ) / (((U-1)+(1-L)*exp(A*x))**2)
  return(list(distance,fprim))
  
}

# TODO : truncated method ?
# TODO : hyperbolic sine

# TODO : shaped -> "param" function which updates parameters is
# contained within distance function
