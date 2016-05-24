require(Matrix)

# PURPOSE: check if a matrix has an intercept in the
#          first row
# ---------------------------------------------------
#  USAGE: has_intercept( x )
#  where  x = n x K matrix 
# ---------------------------------------------------
# RETURN:
#        TRUE - if the matrix has na intercept int the
#               first column
#        FALSE- if the matrix has no intercept 
#        NULL - if the matrix has an intercept in a 
#               column other than the first
has_intercept <- function( x ){
  n         <- nrow(x)
  k         <- ncol(x)
  intercept <- rep(1, n)
  results   <- FALSE
  for( c in 1:k ){
    check <- which( intercept == x[ ,c] )
    if( length( check ) == n && c == 1){
      results <- TRUE
    }else if( length( check ) == n && c != 1){
      print('has_intercept: intercept term must be in first column of the x-matrix')
      results <- NULL
    }
  }
  return( results )
}


sar_bayesian_estimation <- function(results){
  out       <- NULL
  #bayesian estimation
  bout_mean <- as.matrix(c(apply(results$bdraw,2,mean),mean(results$pdraw))) #parameter mean column
  bout_sd   <- as.matrix(c(apply(results$bdraw,2,sd)  ,sd(results$pdraw))) #parameter sd colum
  bout_sig  <- matrix(data=NA, nrow=nrow(bout_mean),ncol=1)
  #build bayesian significance levels
  draws     <- cbind( results$bdraw, results$pdraw )
  for( i in 1:ncol(draws) ){
    if( bout_mean[i,1] > 0){
      cnt <- which( draws[,i] > 0 )
    }else{
      cnt <- which( draws[,i] < 0 )
    }
    bout_sig[i,1] <- 1 - (length(cnt)/(results$ndraw-results$nomit))
  }
  
  out$bhat      <- bout_mean
  out$bhat_sd   <- bout_sd
  out$bhat_Bpval<- bout_sig
  out$bhat_names<- as.matrix( results$names )
  return(out)
}



#----------------------------------------------------
# PURPOSE: compute the eigenvalues for the weight matrix
# ---------------------------------------------------
#  USAGE: results = far_eigs(eflag, W)
#  results.rmin 
#  results.rmax maximum eigen value
#  results.time execution time
# ---------------------------------------------------
# Adapted to R by:
# Miguel Godinho de Matos
# Carnegie Mellon University
# Dept. Engineering & Public Policy
sar_eigs <- function(eflag, W){
  results <- NULL
  results$rmin <- -1
  results$rmax <-  1
  results$time <-  0
  # compute eigenvalues  
  if( eflag == 1 ){
    t0           <- Sys.time()
    lambda       <- as.double( eigen( W, only.values=TRUE )$values )
    results$rmin <- 1/min(lambda)
    results$rmax <- 1 
    results$time <- Sys.time() - t0 
  }
  return( results )
}

#---------------------------------------------------
# PURPOSE: compute the log determinant |I_n - rho*W|
# using the user-selected (or default) method
# ---------------------------------------------------
# USAGE: detval = sar_lndet(lflag,W,rmin,rmax)
# where ldetflag,rmin,rmax,W contains input flags 
# ---------------------------------------------------
# Code is now based on method spdep::do_ldet written by Roger Bivand
# ldetflag = 0 : Pace and Barry (1997) grid
# Pace, R. and Barry, R. (1997). Quick computation of spatial autoregressive estimators.
# Geographics Analysis, 29(3):232-247.
#
# ldetflag = 1 : Pace and LeSage (2004) Chebyshev approximation
# Pace, R. and LeSage, J. (2004). Chebyshev approximation of log-determinants of
# spatial weight matrices. Computational Statistics & Data Analysis, 45(2):179-
# 196.
#
# ldetflag = 2 : Barry and Pace (1999) MC approximation
# Barry, R. and Pace, R. (1999). Monte Carlo estimates of the log determinant of
# large sparse matrices. Linear Algebra and its Applications, 289(1-3):41-54.
sar_lndet <- function(ldetflag,W,rmin,rmax){
  results       <- NULL
  env <- new.env(parent=globalenv())
  listw <- mat2listw(W)
  assign("n", nrow(W), envir=env)
  assign("listw", listw, envir=env)
  assign("family", "SAR", envir=env)
  assign("similar", FALSE, envir=env)

  # do lndet approximation calculations if needed
  if( ldetflag == 0 ){ # no approximation
    t0            <- Sys.time()
    SE_classic_setup(env)
    # fine grid is already computed
    results$detval  <- get("detval1", env)
  }else if( ldetflag == 1 ) { 
    t0   <- Sys.time()
    cheb_setup(env, q=4)  # Chebychev approximation, q = 4
    # For Chebyshev-Approximation we must compute the fine grid with do_ldet
    # to be later used for numerical integration
    detval1 <- seq(rmin, rmax, 0.001)
    detval2 <- sapply(detval1, do_ldet, env)
    results$detval  <- cbind(detval1, detval2)
  } else if( ldetflag == 2 ) {
    t0   <- Sys.time()
    mcdet_setup(env, p=16, m=30, which=1)
    detval1 <- seq(rmin, rmax, 0.001)
    detval2 <- sapply(detval1, do_ldet, env)
    results$detval  <- cbind(detval1, detval2)
  } else{
    # TODO: Pace and Barry, 1998 spline interpolation
    stop('sar_lndet: method not implemented')
  }
  results$time   <- Sys.time() - t0
  return( results )
 }
 


# PURPOSE: computes Pace and Barry's grid for log det(I-rho*W) using sparse matrices
# -----------------------------------------------------------------------
# USAGE: out = lndetfull(W,lmin,lmax)
# where:    
#             lmin  = lower bound on rho
#             lmax  = upper bound on rho
# -----------------------------------------------------------------------
# RETURNS: out = a structure variable
#          out.lndet = a vector of log determinants for 0 < rho < 1
#          out.rho   = a vector of rho values associated with lndet values
# -----------------------------------------------------------------------
# NOTES: should use 1/lambda(max) to 1/lambda(min) for all possible rho values
# -----------------------------------------------------------------------
# References: % R. Kelley Pace and  Ronald Barry. 1997. ``Quick
# Computation of Spatial Autoregressive Estimators'', Geographical Analysis
# -----------------------------------------------------------------------
#Written by:
#James P. LeSage, Dept of Economics
# University of Toledo
# 2801 W. Bancroft St,
# Toledo, OH 43606
# jpl@jpl.econ.utoledo.edu
#
# Adapted to R by:
# Miguel Godinho de Matos
# Carnegie Mellon University
# Dept. Engineering & Public Policy

# @param W spatial weight matrix
# @param rmin
# @param rmax
lndetfull <- function( W, rmin, rmax ){
  rvec    <- seq(rmin, rmax, 0.01)
  In      <- speye(nrow(W))
  niter   <- length(rvec)
  
  results       <- NULL
  results$rho   <- NULL
  results$lndet <- NULL
  
  for(i in 1:niter){
    rho             <- rvec[i]
    z               <- In - rho*W
    results$rho[i]  <- rho
    results$lndet[i]<- as.double(determinant(z, logarithm=TRUE)$modulus) 
  }  
  return(results)
}


# PURPOSE: This function computes the a
# the chebyshev a log determinant approximation. 
# ---------------------------------------------------
#  USAGE: 
#  where:
# ---------------------------------------------------
# SOURCE: 
# Pace and LeSage 2004, Chebyshev Approximation of 
# log-determinants of spatial weight matrices, 
# Computational Statistics and Data Analysis, vol 45
#
# Written by:
# Miguel Godinho de Matos
# Carnegie Mellon University
# Dept. Engineering & Public Policy

lndetChebyshev <- function(W, rmin=-1, rmax=1){
  results       <- NULL
  results$rho   <- NULL
  results$lndet <- NULL
  
  n  <- nrow(W)
  D1 <- W           #keep the notation of the paper
  D2 <- D1 %*% D1
  
  #pre-compute all required traces
  tr <- c(
    n,               #tr(I)
    0,               #tr(D)
    sum(D1 * t(D1)), #tr(D^2)
    sum(D1 * t(D2)), #tr(D^3)
    sum(D2 * t(D2))  #tr(D^4)
  )  
  
  #pre-compute chebyshev approximation
  trTD <- c(
    tr[1],              #tr( TD0 ) = tr(I)
    tr[2],              #tr( TD1 ) = 0
    2*tr[3] - tr[1]  ,  #tr( TD2 ) = 2 * tr(TD2) - tr(I)
    4*tr[4] - 3*tr[2],  #tr( TD3 ) = 4 * tr(TD3) - 3 * tr(TD)
    8*tr[5] - 8*tr[3] + tr[ 1 ]
  ) 
  
  q    <- 4
  rvec <- seq(rmin, rmax, 0.01)
  niter<- length(rvec)
  for( i in 1:niter ){
    alpha <- rvec[ i ]
    c0    <- 0
    for( j in 1:(q+1) ){
      c0 <- c0 + chebyshev(j, q, alpha)*trTD[j]
    }
    results$rho[i]   <- alpha
    results$lndet[i] <- c0 - n/2 * chebyshev(1, q, alpha)
  }
  return( results )
}


# PURPOSE: This function is part of the computation of
# the chebyshev a log determinant approximation. 
# ---------------------------------------------------
#  USAGE: 
#  where:
# ---------------------------------------------------
# SOURCE: 
# Pace and LeSage 2004, Chebyshev Approximation of 
# log-determinants of spatial weight matrices, 
# Computational Statistics and Data Analysis, vol 45
#----------------------------------------------------
# AUTHOR: 
# Miguel Godinho de Matos
# Miguel.GodinhoMatos@gmail.com
# PhD Student of Engineering & Public Policy
# Carnegie Mellon University
chebyshev <- function(j, q, alpha){
  c0 <- 0
  for( k in 1:(q+1)){
    c0 <- c0 + log(1-(alpha*cos((pi*(k-0.5))/(q+1))))*cos((pi*(j-1)*(k-0.5))/(q+1))
  }
  c1 <- (2/(q+1)) * c0
  return(c1)
}
    
# PURPOSE: returns a vector of the log-marginal over a grid of rho-values
# -------------------------------------------------------------------------
# USAGE: out = sar_marginal(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2)
# where:       detval = an ngrid x 2 matrix with rho-values and lndet values
#                  e0 = y - x*b0;
#                 ed = Wy - x*bd;
#               epe0 = e0'*e0;
#               eped = ed'*ed;
#              epe0d = ed'*e0;
#               nobs = # of observations
#               nvar = # of explanatory variables
#            logdetx = log(det(x'*x))
#                 a1 = parameter for beta prior on rho
#                 a2 = parameter for beta prior on rho
# -------------------------------------------------------------------------
# RETURNS: out = a structure variable
#          out = log marginal, a vector the length of detval
# -------------------------------------------------------------------------
# NOTES: -this does not take any prior on beta, sigma into account, 
#        uses diffuse priors 
# -------------------------------------------------------------------------
# written by:
# James P. LeSage, last updated 3/2010
# Dept of Finance & Economics
# Texas State University-San Marcos
# 601 University Drive
# San Marcos, TX 78666
# jlesage@spatial-econometrics.com
#
#Adapted to R by:
# Miguel Godinho de Matos
# Carnegie Mellon University
# Dept. Engineering & Public Policy
sar_marginal <- function(detval,e0,ed,epe0,eped,epe0d,nobs,nvar,logdetx,a1,a2){
  n      <- nrow(detval)
  nmk    <- (nobs - nvar)/2
  bprior <- beta_prior( detval[,1] , a1, a2)
  C      <- log(bprior) + lgamma(nmk) - nmk*log(2*pi) - 0.5 * logdetx
  iota   <- matrix(nrow=n, ncol=1, data=1)
  z      <- epe0 * iota - 2*detval[,1] * epe0d + (detval[,1] * detval[,1]) * eped 
  den    <- C + detval[,2] - nmk*log(z)
  return( base::Re( den ) )
}


# % PURPOSE: computes and prints posterior model probabilities using log-marginals
# % ---------------------------------------------------
# %  USAGE: probs = model_probs(results1,results2,results3, ...)
# %  where: results_matrix is a  matrix of results structures returned by estimation
# %         functions, sar_g, sdm_g, sac_g, sem_g
# % e.g. result1 = sar_g(y,x,W,ndraw,nomit,prior);
# %      result2 = sem_g(y,x,W,ndraw,nomit,prior);
# % results_matrix = [result1 result2];
# % model_probs(results_matrix);
# % ---------------------------------------------------
# %  RETURNS: probs = a vector of posterior model probabilities
# % ---------------------------------------------------
# 
# % written by:
# % James P. LeSage, 7/2003
# % Dept of Economics
# % University of Toledo
# % 2801 W. Bancroft St,
# % Toledo, OH 43606
# % jlesage@spatial-econometrics.com
# Ported to R by Miguel Godinho de Matos
#
# @param detval
# @param result_list a list of estimated model objects
# @return a vector of posterior model probabilities
model_probs <- function(detval, result_list){
  stopifnot(is.list(result_list))
  nmodels   <- length(result_list)
  nrho      <- nrow(detval)
  lmarginal <- matrix(ncol=nmodels, nrow=nrho, data=0)
  
  for(i in 1:nmodels){
    lmarginal[,i] <- result_list[[i]]$mlike
  } 
  
  xx  <- exp(lmarginal- max(lmarginal))
  #yy  <- repmat(as.matrix(detval[,1]), nmodels)
  yy  <- matrix(detval[,1], length(detval[,1]), nmodels)

  isum<- colSums((yy[2:nrho,] + yy[1:(nrho-1),])*(xx[2:nrho,]-xx[1:(nrho-1),])/2)
  psum<- sum( isum )
  return( as.matrix( isum/psum ) )
}

