## Copyright (C) 2005, 2006, 2007, 2008/2006, 2010  Antonio, Fabio Di Narzo
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## A copy of the GNU General Public License is available via WWW at
## http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
## writing to the Free Software Foundation, Inc., 59 Temple Place,
## Suite 330, Boston, MA  02111-1307  USA.

## This program is a translation of Prof. Marcelo Medeiros's Matlab codes
##    and is indebted to him.

calculateLinearCoefficients <- function(A, b){
  #x <- lm.fit(x = A, y = b)$coefficients
  #x <- lm(b ~ . -1, data = data.frame(A))$coefficients  
  x <- as.vector(ginv(A) %*% b)
  x
}

#Logistic transition function
# z: variable
# gamma: smoothing parameter
# th: threshold value
G <- function(z, gamma, th) {
  if((length(th) > 1) && (length(gamma) > 1))
    t( apply( as.matrix(z) , 1, plogis, th, 1/gamma) )
  else
    plogis(z, th, 1/gamma)
}

#Fitted values, given parameters
# phi1: vector of linear parameters
# phi2: vector of tr. functions' parameters
F <- function(phi1, phi2, x_t, s_t) {
  noRegimes <- dim(phi1)[1];

  local <- array(0, c(noRegimes, dim(x_t)[1]))

  local[1,] <- x_t %*% phi1[1,];  
  for (i in 2:noRegimes) 
    local[i,] <-
      (x_t %*% phi1[i,]) * G(s_t, gamma= phi2[i - 1,1], th= phi2[i - 1,2]);

  return(apply(local, 2, sum));
}



#'computeGradient
#'
#'Computes the gradient under the null hypothesis
#'
#'
#'@param object fitted star model
#'@param ... currently unused
#'@return computed gradient
#'@author J. L. Aznarte
#'@seealso \code{\link{addRegime}}
#'@references TODO
#'@keywords ts internal
#'@export
#'@examples
#'
#'##TODO
#'
computeGradient <- function(object, ...)
  UseMethod("computeGradient")

# Computes the gradient
#
# object: a valid STAR model.
#
# Returns a list of the gradients with respect to  the linear and
#     nonlinear parameters
computeGradient.star <- function(object, ...)
{
 
  noRegimes <- object$model.specific$noRegimes;
  gamma <- object$model.specific$phi2[,1];
  th <- object$model.specific$phi2[,2];

  phi1 <- object$model.specific$phi1;

  n.used <- NROW(object$str$xx);
  s_t<- object$model.specific$thVar;
  x_t <- cbind(1, object$str$xx);

  fX <- array(0, c(noRegimes - 1, n.used));
  dfX <- array(0, c(noRegimes - 1, n.used));
  gPhi <- x_t;
  for (i in 1:(noRegimes - 1)) {
    fX[i,] <- G(s_t, gamma[i], th[i]);
    dfX[i,] <- G(s_t, gamma[i], th[i]) * (1 - G(s_t, gamma[i], th[i]));
    gPhi <- cbind(gPhi, kronecker(matrix(1, 1, NCOL(x_t)), fX[i,]) * x_t)
  }
  
  gGamma <- array(0, c(n.used, noRegimes-1));
  gTh <- array(0, c(n.used, noRegimes-1))
  for (i in 1:(noRegimes - 1)) {
    gGamma[, i] <- (x_t %*% phi1[i + 1,]) * (dfX[i,] * (s_t - th[i]));
    gTh[,i] <-         - (x_t %*% phi1[i + 1,]) * (gamma[i] * dfX[i,]);
  }
  
  return(cbind(gPhi, gGamma, gTh))

}

testRegime <- function(object, ...)
  UseMethod("testRegime")

# Tests (within the LM framework), being the null hypothesis H_0 that the
#    model 'object' is complex enough and the alternative H_1 that an extra
#    regime should be considered.
#
# object: a STAR model already built with at least 2 regimes.
# G: the gradient matrix obtained by using computeGradient()
# rob: boolean indicating if robust tests should be used or not.
# sig: significance level of the tests
# 
# returns a list containing the p-value of the F statistic and a boolean,
#      true if there is some remaining nonlinearity and false otherwise.
testRegime.star <- function(object, G, rob=FALSE, sig=0.05, trace = TRUE, ...)
{

  e <-  object$residuals;
  s_t <- object$model.specific$thVar;
  nG <- NCOL(G);
  T <- length(e);

  normG <- norm(Matrix(t(G)) %*% e)
  norm2G <-  sum(abs(t(G) %*% e)^2)^(1/2) 
#  cat("norm2G = ", norm2G, "normG = ", normG, "\n");

  # Compute the rank of G' * G
  G2 <- t(G) %*% G;
  s <- svd(G2);
  tol <- max(dim(G2)) * s[1]$d * 2.2204e-16
  rG <- qr(G2, tol)$rank
#  cat("rango G = ", rG, "\n");
  
#  rG <- sum(svd(G2)$d > tol);

  if (normG > 1e-6) {
    if(rG < nG) {
#      cat("A1\n")
      PCtmp <- princomp(G);
      PC <- PCtmp$loadings;
#      GPCA <- PCtmp$scores;
      lambda <- cumsum(PCtmp$sdev^2 / sum(PCtmp$sdev^2));
      indmin <- min(which(lambda > 0.99999));
      GPCA <- t(PC%*%t(G));
      GPCA <- GPCA[, 1:indmin];
#      b <- solve(t(GPCA) %*% GPCA) %*% t(GPCA) %*% e;
      b <- lm(e ~ . -1, data = data.frame(GPCA))$coefficients
      dim(b) <- c(NCOL(GPCA), 1);
      u <-  e - GPCA %*% b;
      xH0 <- GPCA;
    } else {
#      cat("A2\n")
#      b <- solve(t(G) %*% G) %*% t(G) %*% e;
      b <- lm(e ~ . -1, data=data.frame(G))$coefficients;
      u <- e - G %*% b;
      xH0 <- G;
    }
  } else {
    u <- e;
    if(rG < nG) {
#     cat("B1\n")
      PCtmp <- princomp(G);
      PC <- PCtmp$loadings;
#      GPCA <- PCtmp$scores;
      lambda <- cumsum(PCtmp$sdev^2 / sum(PCtmp$sdev^2));
      indmin <- min(which(lambda> 0.99999));
      GPCA <- t(PC%*%t(G));
      GPCA <- GPCA[, 1:indmin];
      xH0 <- GPCA;
    } else {
#      cat("B2\n")
      xH0 = G;
    }
  }

  SSE0 <- sum(u^2);

  # Regressors under the alternative:
  xx <- cbind(1, object$str$xx);
  nX <- NCOL(xx);
  if (object$model.specific$externThVar) {
    s <- kronecker(matrix(1, 1, nX), s_t); #repmat(s_t, 1, nX)
    xH1 <- cbind(xx * s, xx * (s^2), xx * (s^3))
  } else {
    s <- kronecker(matrix(1, 1, nX - 1), s_t);
    xH1 <- cbind(xx[,2:nX] * s, xx[,2:nX] * (s^2), xx[,2:nX] * (s^3))
  }

  Z <- cbind(xH0, xH1)

  # Standarize the regressors
  nZ <- NCOL(Z);
  sdZ <- apply(Z,2,sd)
  dim(sdZ) <- c(1, nZ)
  sdZ <- kronecker(matrix(1, T, 1), sdZ) # repeat sdZ T rows
  Z[,2:nZ] <- Z[,2:nZ] / sdZ[,2:nZ]

   # Compute the rank of Z
  s <- svd(Z);
  tol <- max(dim(Z)) * s[1]$d * 2.2204e-16
  rZ <- qr(Z, tol)$rank
  if(rZ < NCOL(Z)) warning("Multicollinearity problem. Aborting.\n")

 # Nonlinear model (Alternative hypothesis)
#  c <- solve(t(Z) %*% Z) %*% t(Z) %*% u;
  c <- lm(u ~ . - 1, data=data.frame(Z))$coefficients
  dim(c) <- c(NCOL(Z), 1);
  v <- u - Z %*% c;
  SSE <- sum(v^2);

  # Compute the third order statistic
  nxH0 <- NCOL(xH0);
  nxH1 <- NCOL(xH1);
  
  F = ((SSE0 - SSE) / nxH1) / (SSE / (T - nxH0 - nxH1));

  pValue <- pf(F, nxH1, T - nxH0 - nxH1, lower.tail = FALSE);

  if (pValue >= sig) {
    return(list(remainingNonLinearity = FALSE, pValue = pValue));
  }
  else {
    return(list(remainingNonLinearity = TRUE, pValue = pValue));
  }
}
  


#'addRegime test
#'
#'addRegime test
#'
#'
#'@param object fitted model object with at least 2 regimes
#'@param ... arguments to and from other methods
#'@return A list containing the p-value of the F statistic and a boolean, true
#'if there is some remaining nonlinearity and false otherwise.
#'@author J. L. Aznarte
#'@seealso \code{\link{star}}
#'@references TODO
#'@keywords ts
#'@export
#'@examples
#'
#'##TODO
#'
addRegime <- function(object, ...)
  UseMethod("addRegime")

#' @S3method addRegime star
addRegime.star <- function(object, ...)
{

  noRegimes <- object$model.specific$noRegimes;
  
  gamma <- object$model.specific$phi2[,1];
  th <- object$model.specific$phi2[,2];
  gamma[noRegimes] <- 0.0;
  th[noRegimes] <- 0.0;

  object$model.specific$noRegimes <- object$model.specific$noRegimes + 1;
  object$noRegimes <- object$noRegimes + 1;
  
  object$model.specific$phi2 <- cbind(gamma, th);

  object$model.specific$phi1 <- rbind(object$model.specific$phi1, rnorm(NCOL(object$str$xx) + 1));

  object$model.specific$coefficients <-
                                          c(object$model.specific$phi1, cbind(gamma, th));
  object$coefficients <- c(object$model.specific$phi1, cbind(gamma, th));
  object$model.specific$k <- length(object$coefficients)
  object$k <- length(object$coefficients);

  return(object);
  
}

startingValues <- function(object, ...)
  UseMethod("startingValues")

# Find promising initial values for next regime
#
# object: a valid (estimated) STAR model with m regimes
#
# returns a modified copy of 'object' with good starting values for
#      gamma[noRegime-1] and th[noRegime-1]. It also modifies the linear
#      parameters phi1 (as they are estimated for the new gamma and th).
startingValues.star <- function(object, trace=TRUE, ...)
{

  noRegimes <- object$model.specific$noRegimes;
  gamma <- object$model.specific$phi2[,1];
  th <- object$model.specific$phi2[,2];

  s_t <- object$model.specific$thVar;
  xx <- object$str$xx;
  x_t <- cbind(1, xx)
  yy <- object$str$yy;
  n.used <- NROW(object$str$xx);
  
  bestCost <- Inf;

  # Maximum and minimum values for gamma
  maxGamma <- 40;
  minGamma <- 10;
  rateGamma <- 5; ############################################!!!
  
  # Maximum and minimum values for c
  minTh <- quantile(as.ts(s_t), .1) # percentil 10 of s_t
  maxTh <- quantile(as.ts(s_t), .9) # percentil 90 of s_t
  rateTh <- (maxTh - minTh) / 200; ################################!!!
  
  for(newGamma in seq(minGamma, maxGamma, rateGamma)) {
    for(newTh in seq(minTh, maxTh, rateTh)) {

      gamma[noRegimes - 1] <- newGamma;
      th[noRegimes - 1] <- newTh;

      tmp <- cbind(x_t, matrix(apply(G(s_t, gamma, th), 2, "*",x_t), 
                               nrow = n.used, ncol = (noRegimes - 1) * NCOL(x_t)))
      #newPhi1 <- lm(yy ~ . - 1, data.frame(tmp))$coefficients;
      newPhi1 <- calculateLinearCoefficients(tmp, yy)

      dim(newPhi1) <- c(NCOL(x_t), noRegimes);
      newPhi1 <- t(newPhi1);

      y.hat <- F(newPhi1, cbind(gamma, th), x_t, s_t);
      cost <- crossprod(yy - y.hat)
      
      if(cost <= bestCost) {
        bestCost <- cost;
        bestGamma <- gamma[noRegimes - 1];
        bestTh <- th[noRegimes - 1];
        phi1 <- newPhi1;
      }
    }
  }

  # Reorder the regimes according to the values of th
  if (noRegimes > 2) {
    cat("Reordering regimes...\n")
    th <- sort(th, index.return=TRUE)$x
    ordering <-  sort(th, index.return=TRUE)$ix
    gamma <- gamma[ordering]

    # reestimate phi's
    tmp <- cbind(x_t, matrix(apply(G(s_t, gamma, th), 2, "*",x_t), 
                             nrow = n.used, ncol = (noRegimes - 1) * NCOL(x_t)))
    #phi1<- lm(yy ~ . - 1, as.data.frame(tmp))$coefficients
    phi1 <- calculateLinearCoefficients(tmp, yy)

    dim(phi1) <- c(NCOL(xx) + 1, noRegimes)
    phi1 <- t(phi1)
  }
  else {
    gamma[1] <- bestGamma;
    th[1] <- bestTh;
  }
  
#  if (trace) cat("Starting values fixed for regime ", noRegimes,
#                 "\n   gamma = ", gamma[noRegimes - 1],
#                 ", th = ", th[noRegimes - 1], 
#                 "; SSE = ", bestCost, "\n");

#  if(trace) {
#    cat("   Starting linear values: \n");
#    for(i in 1:noRegimes)
#      cat("   ", phi1[i,], "\n")
#  }


  object$model.specific$phi2 <- cbind(gamma, th);
    
  object$model.specific$phi1 <- phi1;
  object$model.specific$coefficients <- c(phi1, cbind(gamma, th));
  object$coefficients <- c(phi1, cbind(gamma, th));
  object$model.specific$k <- length(object$coefficients)
  object$k <- length(object$coefficients)
  
  return(object);
  
}

estimateParams <- function(object, ...)
  UseMethod("estimateParams")

# Estimates the parameters of a given STAR model.
#
# object: a valid STAR model.
#
# Estimates object$model.specific$phi1
#                   object$model.specific$phi2
estimateParams.star <- function(object, trace=TRUE, control=list(), ...)
{

  s_t<- object$model.specific$thVar;
  xx <- object$str$xx
  x_t <- cbind(1, xx)
  yy <- object$str$yy
  n.used <- NROW(object$str$xx);
  noRegimes <- object$model.specific$noRegimes;

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 2.- Estimate nonlinear parameters
  
  # Function to compute the gradient 
  #
  # Returns the gradient with respect to the error
  gradEhat <- function(phi2, phi1) {
    dim(phi2) <- c(noRegimes - 1, 2)
    gamma <- phi2[,1];
    th <- phi2[,2];

    y.hat <- F(phi1, phi2, x_t, s_t)
    e.hat <- yy - y.hat;

    fX <- array(0, c(noRegimes - 1, n.used));
    dfX <- array(0, c(noRegimes - 1, n.used));
    gPhi <- x_t;
    for (i in 1:(noRegimes - 1)) {
      fX[i,]   <- G(s_t, gamma[i], th[i]);
      dfX[i,] <- G(s_t, gamma[i], th[i]) * (1 - G(s_t, gamma[i], th[i]));
      gPhi <- cbind(gPhi, kronecker(matrix(1, 1, NCOL(x_t)), fX[i,]) * x_t)
    }
    
    gGamma <- array(0, c(n.used, noRegimes-1));
    gTh <- array(0, c(n.used, noRegimes-1))
    for (i in 1:(noRegimes - 1)) {
      gGamma[, i] <- as.vector(x_t %*% phi1[i + 1,]) * as.vector(dfX[i,] *
                                                                 (s_t - th[i]));
      gTh[,i] <-         - as.vector(x_t %*% phi1[i + 1,]) * as.vector(gamma[i] *
                                                                       dfX[i,]);
    }
    
    J = - cbind(gGamma, gTh) / sqrt(n.used)
      
    return(2 * t(e.hat) %*% J)
    
  }

  # Sum of squares function
  # phi2: vector of parameters
  SS <- function(phi2, phi1) {
    dim(phi2) <- c(noRegimes - 1, 2)

    for (i in 1:(noRegimes - 1))
      if(phi2[i,1] <0)
        phi2[i,1] <- - phi2[i,1];
    
    tmp <- cbind(x_t,
                 matrix(apply(G(s_t, gamma = phi2[,1], th = phi2[,2]), 2, "*", x_t),
                             nrow = n.used, ncol = (noRegimes - 1) * NCOL(x_t)))

    #newPhi1<- lm(yy ~ . - 1, as.data.frame(tmp))$coefficients
    newPhi1 <- calculateLinearCoefficients(tmp, yy)

    dim(newPhi1) <- c(NCOL(xx) + 1, noRegimes)
    newPhi1 <- t(newPhi1);

    # Return the sum of squares
    y.hat <- F(newPhi1, phi2, x_t, s_t)
    crossprod(yy - y.hat)
  }
    
  res <- optim(object$model.specific$phi2, SS, #gradEhat,
               control = control, phi1 = object$model.specific$phi1)

  newPhi2 <- res$par
  dim(newPhi2) <- c(noRegimes - 1, 2)
  
  gamma <- newPhi2[,1]
  th <- newPhi2[,2]
  
  for (i in 1:(noRegimes - 1))
    if(gamma[i] <0)
      gamma[i] <- - gamma[i];
  
  if (trace) cat("Optimized values fixed for regime ", noRegimes,
                 ": gamma = ", gamma[noRegimes - 1],
                 ", th = ", th[noRegimes - 1],"\n");

  if(trace) 
    if(res$convergence != 0)
      cat("*** Convergence problem. Code: ",res$convergence,"\n")
    else {
      cat("Optimization algorithm converged\n")
    }

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 3.- Estimate linear parameters again
  tmp <- cbind(x_t, matrix(apply(G(s_t, gamma = gamma, th = th), 2, "*", x_t),
                           nrow = n.used, ncol = (noRegimes - 1) * NCOL(x_t)))
  
  #newPhi1<- lm(yy ~ . - 1, as.data.frame(tmp))$coefficients
  newPhi1 <- calculateLinearCoefficients(tmp, yy)

  dim(newPhi1) <- c(NCOL(xx) + 1, noRegimes)
  newPhi1 <- t(newPhi1)

  if(trace) {
    cat("Optimized linear values: \n");
    for(i in 1:noRegimes)
      cat(newPhi1[i,], "\n")
  }

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 4.- Compute yhat and ehat
  
  object$fitted.values <- F(newPhi1, cbind(gamma, th), x_t, s_t); # y.hat
  object$residuals <- yy - object$fitted.values; # e.hat
  object$coefficients <- c(newPhi1, newPhi2)
  object$model.specific$phi1 <- newPhi1
  object$model.specific$phi2 <- cbind(gamma, th);
  
  return(object)
  
}

# Incremental STAR fitter
#
#   Builds a STAR model with as many regimes as needed, using the
#     bottom-up strategy proposed by Terasvirta et al.
#
#   x: the time series 
#   m: the order of the autoregressive terms
#   d: 
#   steps
#   series
#   rob
#   sig
#' @export
star <- function(x, m=2, noRegimes, d = 1, steps = d, series, rob = FALSE,
                 mTh, thDelay, thVar, sig=0.05, trace=TRUE, control=list(), ...)
{

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 1. Build the nlar object and associated variables.
  if(missing(series))   series <- deparse(substitute(x))

  str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)

  xx <- str$xx
  yy <- str$yy
  externThVar <- FALSE
  n.used <- NROW(xx)

  if(missing(noRegimes))
     noRegimes <- Inf;
    
  if(!missing(thDelay)) {
    
    if(thDelay>=m) 
      stop(paste("thDelay too high: should be < m (=",m,")"))
    
    z <- xx[,thDelay+1]

  } else if(!missing(mTh)) {

    if(length(mTh) != m) 
      stop("length of 'mTh' should be equal to 'm'")

    z <- xx %*% mTh #threshold variable
    dim(x) <- NULL

  }
  else if(!missing(thVar)) {

    if(length(thVar)>nrow(xx)) {
      z <- thVar[1:nrow(xx)]

      if(trace) cat("Using only first", nrow(xx), "elements of thVar\n")

    }
    else
      z <- thVar
    externThVar <- TRUE
  } else {
    if(trace) cat("Using default threshold variable: thDelay=0\n")
    z <- xx[,1]
  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # 2. Linearity testing
  if (trace) cat("Testing linearity...   ")
  testResults <- linearityTest.star(str, z, rob=rob, sig=sig, trace = trace)
  pValue <- testResults$pValue;
  increase <- ! testResults$isLinear;

  if(trace) cat("p-Value = ", pValue,"\n")
  if(!increase) {
    if(trace) cat("The series is linear. Using linear model instead.\n")
    return(linear(x, m=m));
  }
  else {
    if(trace) cat("The series is nonlinear. Incremental building procedure:\n")

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # 3. Build the 2-regime star
    if(trace) cat("Building a 2 regime STAR.\n")
    object <- star.predefined(x, m=m, noRegimes=2, thVar=thVar,
                              d=d, steps=steps, mTh=mTh, thDelay=thDelay,
                              series=series, trace=trace, control=control)
    if(noRegimes==2) {
      if(trace) cat("Finished building a MRSTAR with 2 regimes\n");
      return(object);
    }
    
    if(trace) cat("\nTesting for addition of regime 3.\n");
    
    if(trace) cat("  Estimating gradient matrix...\n");
    G <- computeGradient(object);

    if(trace) cat("  Done. Computing the test statistic...\n");
    testResults <- testRegime(object, G = G, rob = rob, sig = sig, trace=trace);

    increase <- testResults$remainingNonLinearity;
    if(increase && trace) {
      cat("  Done. Regime 3 is needed (p-Value = ", testResults$pValue,").\n");
    }
    else {
      if(trace)
        cat("  Done. Regime 3 is NOT accepted (p-Value = ",
            testResults$pValue,").\n");
    }

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
    # 4. Add-regime loop
    while (increase && (object$model.specific$noRegimes < noRegimes)) {
      
      if(trace) cat("Adding regime ",
                    object$model.specific$noRegimes + 1,".\n");
      object <- addRegime(object);
      
      if(trace) cat("Fixing good starting values for regime ",
                    object$model.specific$noRegimes,"...\n");
      object <- startingValues(object, control=control, trace=trace);
      
      if(trace) cat('Estimating parameters of regime',
                    object$model.specific$noRegimes, '...\n')
      object <- estimateParams(object, control=control, trace=trace);
      
      if(trace) cat("Ok. \n Testing for addition of regime ",
                    object$model.specific$noRegimes + 1, ".\n");
      
      if(trace) cat("  Estimating gradient matrix...\n");
      G <- computeGradient(object);

      if(trace) cat("  Computing the test statistic...\n");
      testResults <- testRegime(object, G = G, rob = rob, sig = sig, trace=trace);

      increase <- testResults$remainingNonLinearity;
      if(increase && trace) {
        cat("Regime ", object$model.specific$noRegimes + 1,
            " is needed (p-Value = ", testResults$pValue,").\n"); 
      }
      else {
        cat("Regime ", object$model.specific$noRegimes + 1,
            " is NOT accepted (p-Value = ", testResults$pValue,").\n");
      }            
    }
    if(trace) cat("\nFinished building a MRSTAR with ",
                  object$model.specific$noRegimes, " regimes\n");
    return(object);
  }
  
}

# Predefined STAR model
#
#   Builds a STAR with a fixed number of regimes ('noRegimes') and
#     fixed parameters phi1 (linear) and phi2 (nonlinear). If no
#     parameters are given they are set to random and evenly
#     distributed, respectively.
#                            TO-DO: proper initial values (grid serch)!!
#
#   mTh, thDelay, thVar: as in setar
#   maxRegressors[i]: maximum number of autoregressors in regime i
#   noRegimes: number of regimes of the model
#   phi1[i]: vector with the maxRegressors[i]+1 parameters of regime i
#   phi2[i]: vector with the parameters of tr. function i.
#   trace: should infos be printed?
#   control: 'control' options to be passed to optim
#
star.predefined <- function(x, m, noRegimes, d=1, steps=d, series,
                  maxRegressors, phi1, phi2,
                  mTh, thDelay, thVar, trace=TRUE, control=list())
{

  if(noRegimes == 1) 
   stop("A STAR with 1 regime is an AR model: use the linear model instead.")
  
  if(noRegimes == 2)
  {

    temp <- lstar(x, m=m, d=d, steps=steps, series=series, mTh=mTh,
                  mH=m, mL=m, thDelay=thDelay, thVar=thVar,
                  trace=trace, control=list())

  ## thVar contains starting NA values: remove:
    temp$model.specific$thVar <-  tail(temp$model.specific$thVar, -(temp$str$n.used-nrow(temp$str$xx)))

    # Cast the lstar into a valid star object
    temp$model.specific$noRegimes <- noRegimes;
    temp$model.specific$m <- m;

#    if(missing(thVar))
#      temp$model.specific$externThVar <- FALSE
#    else temp$model.specific$externThVar <- TRUE

    temp$model.specific$phi1 <- # Linear parameters
      temp$model.specific$coefficients[1:((m + 1) * 2)];
    dim(temp$model.specific$phi1) <- c(m + 1, noRegimes)
    temp$model.specific$phi1 <- t(temp$model.specific$phi1)

    temp$model.specific$phi2 <- # Non linear parameters
      temp$model.specific$coefficients[(((m + 1) * 2) + 1):
                                       (((m + 1) * 2) + 2)];
    dim(temp$model.specific$phi2) <- c(1, 2)

#if(trace) cat("Optimized values: gamma = ",temp$model.specific$phi2[1],
#                  ", th =", temp$model.specific$phi2[2], "\n")
#  if(trace) cat("Optimized linear values:\n", temp$model.specific$phi1[1,],
#                "\n", temp$model.specific$phi1[2,], "\n")
    
    return(extend(nlar(str=temp$str, 
                       coefficients = temp$model.specific$coefficients,
                       fitted.values = temp$fitted,
                       residuals = temp$residuals,
                       noRegimes = noRegimes,
                       k = temp$k,
		       model=temp$model,
                       model.specific = temp$model.specific), "star"))
  } 
  
  if(missing(m))
    m <- max(maxRegressors, thDelay+1)

  if(missing(series))
    series <- deparse(substitute(x))

  str <- nlar.struct(x=x, m=m, d=d, steps=steps, series=series)
  xx <- str$xx
  yy <- str$yy
  externThVar <- FALSE
  n.used <- NROW(xx)

  if (missing(maxRegressors)) {
    maxRegressors <- rep(m, times = noRegimes)
    if (trace) 
      cat("Using maximum autoregressive order for all regimes: ", m,"\n")
  }

  if(!missing(thDelay))
  {
    if(thDelay>=m) 
      stop(paste("thDelay too high: should be < m (=",m,")"))
    z <- xx[,thDelay+1]
  } else if(!missing(mTh))
  {
    if(length(mTh) != m) 
      stop("length of 'mTh' should be equal to 'm'")
    z <- xx %*% mTh #threshold variable
    dim(x) <- NULL
  }
  else if(!missing(thVar))
  {
    if(length(thVar)>nrow(xx)) {
      z <- thVar[1:nrow(xx)]
      if(trace) 
        cat("Using only first", nrow(xx), "elements of thVar\n")
    }
    else 
      z <- thVar
    externThVar <- TRUE
  } else
  {
    if(trace) 
      cat("Using default threshold variable: thDelay=0\n")
    z <- xx[,1]
  }
  
#Automatic starting values####################
# TO DO: grid search over phi2.

  if(missing(phi2)) {         
    phi2 <- array(0, c(noRegimes - 1, 2))
    range <- range(z)
    phi2[,1] <- 0.5                   # phi2[,1] = gamma
    phi2[1:(noRegimes - 1),2] <- # phi2[,2] = th
      range[1] + ((range[2] - range[1]) / (noRegimes - 1)) * (0:(noRegimes-2))
    if(trace) {
      cat('Missing starting transition values. Using uniform distribution.\n')
    }
    optimize <- TRUE
  }
  else {
    optimize <- FALSE
  }

  if(missing(phi1)) {
    phi1 <- array(rnorm(noRegimes * m+1), c(noRegimes, m+1))
    if(trace) {
      cat('Missing starting linear values. Using random values.\n')
    }
  }

#############################################

  #Fitted values, given parameters
  # phi1: vector of linear parameters
  # phi2: vector of tr. functions' parameters
  F <- function(phi1, phi2) {
    x_t <- cbind(1, xx)
    local <- array(0, c(noRegimes, n.used))
    local[1,] <- x_t %*% phi1[1,];

    for (i in 2:noRegimes) 
      local[i,] <-
        (x_t %*% phi1[i,]) * G(z, gamma= phi2[i - 1,1], th= phi2[i - 1,2]);
    
    result <- apply(local, 2, sum)
    result
  }
  
  #Sum of squares function
  #p: vector of parameters
  SS <- function(phi2) {
    dim(phi2) <- c(noRegimes - 1, 2)
    
    # We fix the linear parameters here, before optimizing the nonlinear.
    tmp <- rep(cbind(1,xx), noRegimes)
    dim(tmp) <- c(NROW(xx), NCOL(xx) + 1, noRegimes)
    for (i in 2:noRegimes) 
      tmp[,,i] <- tmp[,,i] * G(z, gamma = phi2[i - 1,1], th = phi2[i - 1,2])

    # Compute the rank of tmp
    s <- svd(tmp);
    tol <- max(dim(tmp)) * s[1]$d * 2.2204e-16
    rtmp <- qr(tmp, tol)$rank
    if(rtmp < NCOL(tmp))
      stop("Multicollinearity problem. Aborting.\n")
      
    #newPhi1<- lm(yy ~ . - 1, as.data.frame(tmp))$coefficients
    newPhi1 <- calculateLinearCoefficients(tmp, yy)

    dim(newPhi1) <- c(noRegimes, m + 1)

    # Return the sum of squares
    y.hat <- F(newPhi1, phi2)
    crossprod(yy - y.hat)
  }
  
#Numerical optimization##########
  res <- list()
  res$convergence <- NA
  res$hessian <- NA
  res$message <- NA
  res$value <- NA

  if(optimize) {
    if(trace) cat('Optimizing...')
    
    p <- as.vector(phi2)
    
    res <- optim(p, SS, hessian = TRUE, control = control)
    
    if(trace) 
      cat(' Nonlinear optimisation done.')

    phi2 <- res$par
    dim(phi2) <- c(noRegimes - 1,2)
    
    for (i in 1:(noRegimes - 1))
      if(phi2[i,1] <0)
        phi2[i,1] <- - phi2[i,1];

#    x_t <- cbind(1,xx);
#    local <- array(0, c(noRegimes, n.used))
#    local[1,] <- x_t %*% phi1[1,]
#    for (i in 2:noRegimes) 
#      local[i,] <-
#        (x_t %*% phi1[i,]) * sigmoid(phi2[i - 1,1] * (z - phi2[i - 1, 2]))
    
#    phi1<- lm(yy ~ . - 1, as.data.frame(local))$coefficients
#    dim(phi1) <- c(noRegimes, m+1)
    
    # We fix the linear parameters here, before optimizing the nonlinear.
    tmp <- rep(cbind(1,xx), noRegimes)
    dim(tmp) <- c(NROW(xx), NCOL(xx) + 1, noRegimes)
    for (i in 2:noRegimes) 
      tmp[,,i] <- tmp[,,i] * G(z, gamma = phi2[i - 1,1], th = phi2[i - 1,2])

    #phi1<- lm(yy ~ . - 1, as.data.frame(tmp))$coefficients
    phi1 <- calculateLinearCoefficients(tmp, yy)

    dim(phi1) <- c(noRegimes, m + 1)

    if(trace) 
      cat(' Linear optimisation done.\n')

    if(trace)
      if(res$convergence!=0)
        cat("Convergence problem. Convergence code: ",res$convergence,"\n")
      else
       cat("Optimization algorithm converged\n")
  }
################################
  
  #Results storing################
  res$phi1 <- phi1
  res$phi2 <- phi2
  
  res$externThVar <- externThVar

  if(!externThVar) {
    if(missing(mTh)) {
      mTh <- rep(0,m)
      mTh[thDelay+1] <- 1
    }
    res$mTh <- mTh
  }
  
  res$thVar <- z
  
  fitted <- F(phi1, phi2)
  
  residuals <- yy - fitted
  dim(residuals) <- NULL	#this should be a vector, not a matrix

  res$noRegimes <- noRegimes
  res$m <- m;
  
  if(!optimize) {
    res$convergence=NA
    res$hessian=NA
    res$message=NA
    res$value=NA
    res$fitted <- fitted
    res$residuals <- residuals
    dim(res$residuals) <- NULL #this should be a vector, not a matrix
    res$k <- length(as.vector(phi1)) + length(as.vector(phi2))
  }
  
  return(extend(nlar(str, 
                     coefficients= c(phi1, phi2),
                     fitted.values = fitted,
                     residuals = residuals,
                     k   = length(as.vector(phi1)) +
                                 length(as.vector(phi2)),
		    model=NULL,
                     model.specific=res), "star"))
}

oneStep.star <- function(object, newdata, itime, thVar, ...)
{

  noRegimes <- object$model.specific$noRegimes;

  phi1 <- object$model.specific$phi1;
  phi2 <- object$model.specific$phi2;

  if(object$model.specific$externThVar) {
    z <- thVar[itime]
  } else {
    z <- newdata %*% object$model.specific$mTh;
    dim(z) <- NULL;
  }

  if(nrow(newdata) > 1) {
    accum <- array(0, c(noRegimes, nrow(newdata)))
    accum[1,] <- (cbind(1,newdata) %*% phi1[1,]);
    for (i in 2:noRegimes) 
      accum[i,] <-
        (cbind(1,newdata) %*% phi1[i,]) * G(z, phi2[i - 1,1], phi2[i - 1,2])
    
    result <- apply(accum, 2, sum)
  }
  else {
    accum <- c(1, newdata) %*% phi1[1,];
    for (i in 2:noRegimes) 
      accum <- accum + 
        (c(1, newdata) %*% phi1[i,]) * G(z, phi2[i - 1,1], phi2[i - 1,2])
    
    result <- accum
  }

  result
  
}

#' @S3method print star
print.star <- function(x, ...) {
  NextMethod(...)
  cat("\nMultiple regime STAR model\n\n")
  x <- x$model.specific
  dg <- options()$digits
  for (i in 1:x$noRegimes) {
    cat("Regime ", i, ":\n")

    cat("    Linear parameters: ")
    cat(paste(round(x$phi1[i,],dg), collapse=", "), '\n')
    if(i > 1) {
      cat("    Non-linear parameters:\n")
      cat(paste(round(x$phi2[i-1,],dg), collapse=", "))
    }
    cat("\n")

  }
  
  invisible(x)
}

