# copyright (C) 2014-2016 A.Rebecq
### Partial calibration functions. Are all private, used by the main
### "calibration" function

penalizedCalib <- function(Xs, d, total, q=NULL, method=NULL, bounds = NULL,
                  alpha = NULL, costs, uCostPenalized=1e2,
                  maxIter=500, calibTolerance=1e-06, lambda=NULL, gap=NULL) {

  ## No distance implemented in penalizedCalib requires updateParameters or params
  updateParameters <- NULL
  params <- NULL
  
  setLambdaPerso <- FALSE
  if(!is.null(lambda)) {
    setLambdaPerso <- TRUE
  }
  
  return(penalCalibAlgorithm(Xs, d, total, q, method,
                        updateParameters, params, costs, uCostPenalized=uCostPenalized
                        , maxIter, calibTolerance, lambda=lambda, setLambdaPerso=setLambdaPerso, gap=gap))

}

penalCalibAlgorithm <- function(Xs, d, total, q=NULL,
                                method, updateParameters, params, costs, uCostPenalized=1e2
                                , maxIter=500, calibTolerance=1e-06
                                , lambda=NULL, setLambdaPerso=FALSE, gap=NULL) {

  ## TODO : change if method is not linear
  distance <- distanceKhiTwo
  
  if(is.null(q)) {
    q <- rep(1,length(d))
  }

  ## In this case, penalized calibration has to check with regular calibration
  checkInfiniteCosts <- cleanCosts(costs, uCostPenalized) == cleanCosts(rep(Inf, length(costs)), uCostPenalized)
  if( all(checkInfiniteCosts) ) {
    matchClassicCalibration <- TRUE
  } else {
    matchClassicCalibration <- FALSE
  }

  if(!is.null(costs)) {
      costs <- cleanCosts(costs, uCostPenalized)
    } else {
      stop("Can't process NULL costs !")
  }

  wRegularCalibration <- NULL
  wRegularCalibration <- calib(Xs, d, total, q, "linear")*d
  
  if(!is.null(wRegularCalibration) && matchClassicCalibration) {
    warning("All costs are infinite: return regular calibration with linear distance")
    return(wRegularCalibration)
  }
  
  ## TODO : handle case when no convergence for linear calib (fail)
  ## (cannot process the following steps)
  
  costs_test <- costs
  costs_test[is.infinite(costs_test)] <- 1e9
  lambdaTest <- distance(wRegularCalibration,d, params) / distanceKhiTwo(costs_test*(d %*% Xs), costs_test*total)

  if(lambdaTest == 0) {
    lambdaTest <- 1
  }

  # If lambda is NULL, lambda is computed
  # so that the two terms of the optimization program
  # are somewhat numerically equal (distance between w and d is
  # estimated by distance on regular linear calibration)
  if(is.null(lambda)) {

    if(matchClassicCalibration) {
#       lambda <- lambdaTest*1e15
      lambda <- 1
    } else {
      lambda <- lambdaTest*1e15
    }

  }

  if(is.na(lambda)) {
    lambda <- 1
  }

  if(!is.null(gap)) {

    ## Check if regular linear calibration weights' gap
    ## is inferior to selected gap
    ## (then it is better to return linear weights)
    gTest <- wRegularCalibration/d

    if( abs(max(gTest) - min(gTest) ) <= gap) {
      warning("Return linear calibration weights, which give inferior gap.")
      return(wRegularCalibration)
    }

    return( searchLambda(Xs, d, total, q,
                                     method, updateParameters, params, costs,
                                     maxIter, calibTolerance
                                      , lastLambdas=NULL, lambdaBackup=1
                                      , lambdaTest=lambda, setLambdaPerso=setLambdaPerso
                                      , gap=gap) )
  }


  paramInit <- wRegularCalibration
  if(is.null(paramInit) || length(paramInit) == 0) {
    paramInit <- d
  }


  if(identical(method,"linear")) {
    ## Solve with linear method (analytical solution)
    A <- t(Xs * d * q) %*% Xs
    C_m_inv <- diag(1/costs)
    w_solution <- d + d*q* (Xs %*% ( ginv(A + lambda*C_m_inv) %*% t(unname(total - d%*%Xs)) ))
  } else {
    
    ## Solve for raking method (entropy distance) by ICRS algorithm
    wTemp <- paramInit
    A <- t(Xs * d * q) %*% Xs
    C_m_inv <- diag(1/costs)
    cont <- TRUE
    while(cont) {
      wTempBackup <- wTemp
      qTemp <- as.vector(2*(wTemp - d) / (d * distanceRaking(wTemp,d)[[2]]))
      
      qTemp[is.na(qTemp)] <- 0
      qTemp[is.infinite(qTemp)] <- 0
      
      A <- t(Xs * d * q*qTemp) %*% Xs
      wTemp <- d + d*q*qTemp* (Xs %*% ( ginv(A + lambda*C_m_inv) %*% t(unname(total - d%*%Xs)) ))
      
      if( all(abs(wTemp - wTempBackup)) <= calibTolerance ) {
        cont <- FALSE
        w_solution <- wTemp
      }
    }
    
  }
  return(w_solution)
}

toOptimize <- function(w, Xs,d, total, lambda, costs, distance, params=NULL) {

  totalsW <- w %*% Xs

  returnD <- distance(w,d, params) + lambda * distanceKhiTwo(costs*totalsW, costs*total)


  return(returnD)
}

cleanCosts <- function(costs, uCostPenalized=1e-2) {
  replacedCosts <- costs
  # replacedCosts[is.infinite(replacedCosts) | replacedCosts < 0] <- infinity
  replacedCosts[replacedCosts < 0] <- Inf

  replacedCosts <- replacedCosts*uCostPenalized

  return(replacedCosts)
}

# @param lastLambdas = c(lastLambdaTooSmall, lastLambdaTooBig)
searchLambda <- function(Xs, d, total, q=rep(1,length(d)),
                         method, updateParameters, params, costs,
                         maxIter=500, calibTolerance=1e-06
                         , lastLambdas=NULL, lambdaBackup=1, lambdaTest=NULL
                         , setLambdaPerso=FALSE
                         , gap, count=0) {
  
  ## TODO : change if method is not linear
  distance <- distanceKhiTwo
  
  # In pratice, optimal lambdas can be found in a very wide domain,
  # so we must set very wide default values
  defaultLastLambdas <- c(1e-50,1e50)

  if(is.null(lambdaTest) && is.null(lastLambdas)) {
    lambdaTest <- 1
    lastLambdas <- defaultLastLambdas
  } else {

    if(is.null(lambdaTest) && !is.null(lastLambdas)) {
      # Don't let search stuck in too great values
      if( abs(lastLambdas[2] - lastLambdas[1]) >= 1e5 ) {
        lambdaTest <- 10**(mean(log(lastLambdas, 10)))
      } else {
        lambdaTest <- mean(lastLambdas)
      }
    }

    if(!is.null(lambdaTest) && is.null(lastLambdas)) {
        lastLambdas <- c(lambdaTest/1e15, lambdaTest*1e15)

        if(setLambdaPerso) {
          warning("Personal lambda set : search interval reduced, may not converge")
          lastLambdas <- c(lambdaTest/1e3, lambdaTest*1e3)
        }
    }

  }

  writeLines(paste("Test with lambda = ", lambdaTest, sep=""))

  w <- penalCalibAlgorithm(Xs, d, total, q,
                                       method, updateParameters, params, costs,
                                       maxIter, calibTolerance
                                       , lambda=lambdaTest, setLambdaPerso=FALSE)

  g <- w/d
  minG <- min(g)
  maxG <- max(g)
  print(minG)
  print(maxG)

  ## Due to the huge numeric range of penalized calibration problems,
  ## we may have to tweak the parameter lambdaTolerance a little before
  ## finding the "right value"
  lambdaTolerance <- 1e-20
  if( (count > 0 && abs(lambdaBackup - lambdaTest) <= lambdaTolerance) || (abs(maxG - minG - gap) <= 1e-4) ) {
    writeLines(paste("Found lambda = ",lambdaTest, " ; count = ",count, sep=""))
    return(w)
  }

  if(maxG - minG >= gap) {
    writeLines("Lambda too small")
    return( searchLambda(Xs, d, total, q=rep(1,length(d)),
                                     method, updateParameters, params, costs,
                                     maxIter=500, calibTolerance=1e-06
                                     , lastLambdas=c(lambdaTest, lastLambdas[2]), lambdaBackup=lambdaTest
                                      , gap=gap, count=count+1) )
  } else {
    writeLines("Lambda too big")
    return( searchLambda(Xs, d, total, q=rep(1,length(d)),
                         method, updateParameters, params, costs,
                         maxIter=500, calibTolerance=1e-06
                         , lastLambdas=c(lastLambdas[1], lambdaTest), lambdaBackup=lambdaTest
                         , gap=gap, count=count+1) )
  }

}


formatCosts <- function(costs, marginMatrix, popTotal) {

  costsFormatted <- NULL

  for(i in 1:nrow(marginMatrix)) {

    nModalitiesMargin <- as.numeric(marginMatrix[i,2])

    if( nModalitiesMargin <= 1 ) {
      costsFormatted <- c(costsFormatted, costs[i])
    } else {
      costsFormatted <- c(costsFormatted, rep(costs[i],nModalitiesMargin))
    }

  }

  ## Add popTotal if not NULL (calibration always exact on this margin)
  if( !is.null(popTotal) ) {
    costsFormatted <- c(costsFormatted, Inf)
  }

  return(costsFormatted)

}

###### Distances #######

# This one is a real distance (between two vectors)
distanceKhiTwo <- function(w,d, params=NULL) {
  return( sum((w-d)*(w-d) / d) )
}

# params are bounds
distanceLogitBounds <- function(w,d, bounds) {

  if(length(bounds) != 2) {
    stop("Must enter LO and UP bounds in a vector.
          Example : bounds=c(0.5,1.5) for LO = 0.5 and UP=1.5")
  }

  L <- bounds[1]
  U <- bounds[2]
  A <- (U - L) / ((1-L)*(U-1))
  r <- w/d
  distanceVec <- 1/A*( (r-L)*log((r-L)/(1-L)) + (U-r)*log((U-r)/(U-1)))
  distance <- d*distanceVec

  return( sum(distance) )
}

distanceTruncated <- function(w,d,bounds, infinity=1e10) {

  r <- w/d
  L <- bounds[1]
  U <- bounds[2]

  r[r > U] <- infinity
  r[r < L] <- infinity

  distanceVec <- 0.5*(r-1)*(r-1)
  distance <- d*distanceVec

  return(sum(distance))
}

## List with distance and its derivative
distanceRaking <- function(w,d, params=NULL) {

  distanceVec <- (w/d)*log(w/d) - (w/d) + 1
  distance <- d*distanceVec
  
  distancePrime <- log(w/d)

  return( list(sum(distance), sum(distancePrime)) )
}

