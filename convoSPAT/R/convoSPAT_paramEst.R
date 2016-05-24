#======================================================================================
# Local likelihood estimation for covariance functions with spatially-varying
# parameters: the convoSPAT() package for R
# Mark D. Risser / The Ohio State University / 2014-2015
#======================================================================================

#======================================================================================
# Local/Global Parameter Estimation
#======================================================================================

#======================================================================================
# Global parameter estimation
#======================================================================================
# For a given correlation matrix, the following functions calculate the global
# parameters based on what is NOT specified to be spatially-varying.
#======================================================================================

#===================================
# First, models without kappa
# Estimates: nugget, variance

#ROxygen comments ----
#' Constructor functions for global parameter estimation.
#'
#' This function generates another function to be used within \code{optim} to
#' obtain maximum likelihood estimates of
#' global variance parameters tausq, sigmasq with a fixed correlation
#' matrix (smoothness is fixed).
#'
#' @param data A vector or matrix of data to use in the likelihood
#' calculation.
#' @param Xmat The design matrix for the mean model.
#' @param Dmat A diagonal matrix of eigenvalues for the covariance
#' matrix.
#' @param Vmat An orthogonal matrix of eigenvectors for the covariance
#' matrix.
#'
#' @return This function returns another function for use in \code{optim}.
#'
#' @examples
#' \dontrun{
#' make_global_loglik1( data, Xmat, Dmat, Vmat )
#' }
#'
#' @export

make_global_loglik1 <- function( data, Xmat, Dmat, Vmat ){

  fixed <- c(FALSE, FALSE)
  params <- fixed
  function(p){
    params[!fixed] <- p

    # Parameters
    tausq <- params[1]
    sigmasq <- params[2]

    N <- dim(data)[1]
    p <- dim(data)[2]

    Cov.inv <- (1/tausq)*( diag(rep(1,N)) - (sigmasq/tausq)*Vmat
                           %*% diag( 1/(1/diag(Dmat) + rep(sigmasq/tausq,N)) )%*%t(Vmat) )

    XCX <- t(Xmat) %*% Cov.inv %*% Xmat
    XCX.inv <- chol2inv( chol(XCX) )

    Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

    # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
    log.det.Cov <- sum( log( 1/diag(Dmat) + rep(sigmasq/tausq,N) ) ) + sum( log(diag(Dmat)) ) + N*log(tausq)
    loglikelihood <- 0.5*( p*log.det.Cov + p*log(det(XCX)) + sum(diag(t(data)%*%Ptemp%*%data)) )

    if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
    return(loglikelihood)

  }
}

# Estimates: variance

#ROxygen comments ----
#' Constructor functions for global parameter estimation.
#'
#' This function generates another function to be used within \code{optim} to
#' obtain maximum likelihood estimates of
#' global variance parameter sigmasq with a fixed correlation
#' matrix (smoothness is fixed). The nugget variance is taken
#' to be spatially-varing.
#'
#' @param data A vector or matrix of data to use in the likelihood
#' calculation.
#' @param Xmat The design matrix for the mean model.
#' @param Corr The correlation matrix matrix.
#' @param obs.nuggets A vector containing the spatially-varying nuggets
#' corresponding to each data location.
#'
#' @return This function returns another function for use in \code{optim}.
#'
#' @examples
#' \dontrun{
#' make_global_loglik2( data, Xmat, Corr, obs.nuggets )
#' }
#'
#' @export

make_global_loglik2 <- function( data, Xmat, Corr, obs.nuggets ){

  fixed <- c(FALSE)
  params <- fixed
  function(p){
    params[!fixed] <- p

    # Parameters
    sigmasq <- params[1]

    N <- dim(data)[1]
    p <- dim(data)[2]

    Edecomp <- eigen(sigmasq*Corr + diag(obs.nuggets))
    Dmat.temp <- diag(Edecomp$values)
    Vmat.temp <- Edecomp$vectors

    Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

    XCX <- t(Xmat) %*% Cov.inv %*% Xmat
    XCX.inv <- chol2inv( chol(XCX) )

    Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

    # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
    loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(data)%*%Ptemp%*%data)) )

    if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
    return(loglikelihood)

  }
}


# Estimates: nugget

#ROxygen comments ----
#' Constructor functions for global parameter estimation.
#'
#' This function generates another function to be used within \code{optim} to
#' obtain maximum likelihood estimates of
#' global variance parameter tausq with a fixed correlation
#' matrix (smoothness is fixed). The process variance is taken
#' to be spatially-varing.
#'
#' @param data A vector or matrix of data to use in the likelihood
#' calculation.
#' @param Xmat The design matrix for the mean model.
#' @param Corr The correlation matrix matrix.
#' @param obs.variance A vector containing the spatially-varying variance
#' corresponding to each data location.
#'
#' @return This function returns another function for use in \code{optim}.
#'
#' @examples
#' \dontrun{
#' make_global_loglik3( data, Xmat, Corr, obs.variance )
#' }
#'
#' @export

make_global_loglik3 <- function( data, Xmat, Corr, obs.variance ){

  fixed <- c(FALSE)
  params <- fixed
  function(p){
    params[!fixed] <- p

    tausq <- params[1]

    N <- dim(data)[1]
    p <- dim(data)[2]

    Edecomp <- eigen(diag(sqrt(obs.variance)) %*% Corr %*% diag(sqrt(obs.variance)) + diag(rep(tausq,N)))
    Dmat.temp <- diag(Edecomp$values)
    Vmat.temp <- Edecomp$vectors

    Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

    XCX <- t(Xmat) %*% Cov.inv %*% Xmat
    XCX.inv <- chol2inv( chol(XCX) )

    Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

    # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
    loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(data)%*%Ptemp%*%data)) )

    if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
    return(loglikelihood)

  }
}

#===================================

#===================================
# Next, models with kappa
# Estimates: nugget, variance, kappa

#ROxygen comments ----
#' Constructor functions for global parameter estimation.
#'
#' This function generates another function to be used within \code{optim} to
#' obtain maximum likelihood estimates of
#' global variance parameters tausq, sigmasq, and nu.
#'
#' @param data A vector or matrix of data to use in the likelihood
#' calculation.
#' @param Xmat The design matrix for the mean model.
#' @param cov.model String; the covariance model.
#' @param Scalemat Matrix; contains the scaling quantities from the
#' covariance function.
#' @param Distmat Matrix; contains the scaled distances.
#'
#' @return This function returns another function for use in \code{optim}.
#'
#' @examples
#' \dontrun{
#' make_global_loglik1_kappa( data, Xmat, cov.model, Scalemat, Distmat )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial

make_global_loglik1_kappa <- function( data, Xmat, cov.model, Scalemat, Distmat ){

  fixed <- c(FALSE, FALSE, FALSE)
  params <- fixed
  function(p){
    params[!fixed] <- p

    # Parameters
    tausq <- params[1]
    sigmasq <- params[2]
    kapp <- params[3]

    Unscl.corr <- cov.spatial( Distmat, cov.model = cov.model,
                               cov.pars = c(1,1), kappa = kapp )
    Corr <- Scalemat*Unscl.corr

    N <- dim(data)[1]
    p <- dim(data)[2]

    Edecomp <- eigen(sigmasq*Corr + diag(rep(tausq,N)))
    Dmat.temp <- diag(Edecomp$values)
    Vmat.temp <- Edecomp$vectors

    Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

    XCX <- t(Xmat) %*% Cov.inv %*% Xmat
    XCX.inv <- chol2inv( chol(XCX) )

    Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

    # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
    loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(data)%*%Ptemp%*%data)) )

    if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
    return(loglikelihood)

  }
}


# Estimates: variance, kappa
#ROxygen comments ----
#' Constructor functions for global parameter estimation.
#'
#' This function generates another function to be used within \code{optim} to
#' obtain maximum likelihood estimates of
#' global variance parameters sigmasq and nu. The nugget variance is
#' taken to be spatially-varying.
#'
#' @param data A vector or matrix of data to use in the likelihood
#' calculation.
#' @param Xmat The design matrix for the mean model.
#' @param cov.model String; the covariance model.
#' @param Scalemat Matrix; contains the scaling quantities from the
#' covariance function.
#' @param Distmat Matrix; contains the scaled distances.
#' @param obs.nuggets A vector containing the spatially-varying nuggets
#' corresponding to each data location.
#'
#' @return This function returns another function for use in \code{optim}.
#'
#' @examples
#' \dontrun{
#' make_global_loglik2_kappa( data, Xmat, cov.model, Scalemat, Distmat, obs.nuggets )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial

make_global_loglik2_kappa <- function( data, Xmat, cov.model, Scalemat, Distmat, obs.nuggets ){

  fixed <- c(FALSE, FALSE)
  params <- fixed
  function(p){
    params[!fixed] <- p

    # Parameters
    sigmasq <- params[1]
    kapp <- params[2]

    Unscl.corr <- cov.spatial( Distmat, cov.model = cov.model,
                               cov.pars = c(1,1), kappa = kapp )
    Corr <- Scalemat*Unscl.corr

    N <- dim(data)[1]
    p <- dim(data)[2]

    Edecomp <- eigen(sigmasq * Corr + diag(obs.nuggets))
    Dmat.temp <- diag(Edecomp$values)
    Vmat.temp <- Edecomp$vectors

    Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

    XCX <- t(Xmat) %*% Cov.inv %*% Xmat
    XCX.inv <- chol2inv( chol(XCX) )

    Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

    # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
    loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(data)%*%Ptemp%*%data)) )

    if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
    return(loglikelihood)

  }
}

# Estimates: nugget, kappa
#ROxygen comments ----
#' Constructor functions for global parameter estimation.
#'
#' This function generates another function to be used within \code{optim} to
#' obtain maximum likelihood estimates of
#' global variance parameters tausq and nu. The process variance is
#' taken to be spatially-varying.
#'
#' @param data A vector or matrix of data to use in the likelihood
#' calculation.
#' @param Xmat The design matrix for the mean model.
#' @param cov.model String; the covariance model.
#' @param Scalemat Matrix; contains the scaling quantities from the
#' covariance function.
#' @param Distmat Matrix; contains the scaled distances.
#' @param obs.variance A vector containing the spatially-varying variance
#' corresponding to each data location.
#'
#' @return This function returns another function for use in \code{optim}.
#'
#' @examples
#' \dontrun{
#' make_global_loglik3_kappa( data, Xmat, cov.model, Scalemat, Distmat, obs.variance )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial

make_global_loglik3_kappa <- function( data, Xmat, cov.model, Scalemat, Distmat, obs.variance ){

  fixed <- c(FALSE, FALSE)
  params <- fixed
  function(p){
    params[!fixed] <- p

    # Parameters
    tausq <- params[1]
    kapp <- params[2]

    Unscl.corr <- cov.spatial( Distmat, cov.model = cov.model,
                               cov.pars = c(1,1), kappa = kapp )
    Corr <- Scalemat*Unscl.corr

    N <- dim(data)[1]
    p <- dim(data)[2]

    Edecomp <- eigen(diag( sqrt(obs.variance) ) %*% Corr %*% diag( sqrt(obs.variance) ) + tausq*diag(rep(1,N)))
    Dmat.temp <- diag(Edecomp$values)
    Vmat.temp <- Edecomp$vectors

    Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

    XCX <- t(Xmat) %*% Cov.inv %*% Xmat
    XCX.inv <- chol2inv( chol(XCX) )

    Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

    # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
    loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(data)%*%Ptemp%*%data)) )

    if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
    return(loglikelihood)

  }
}


# Estimates: kappa
#ROxygen comments ----
#' Constructor functions for global parameter estimation.
#'
#' This function generates another function to be used within \code{optim} to
#' obtain maximum likelihood estimates of
#' global variance parameters nu. The process variance
#' and nugget variance are taken to be spatially-varying.
#'
#' @param data A vector or matrix of data to use in the likelihood
#' calculation.
#' @param Xmat The design matrix for the mean model.
#' @param cov.model String; the covariance model.
#' @param Scalemat Matrix; contains the scaling quantities from the
#' covariance function.
#' @param Distmat Matrix; contains the scaled distances.
#' @param obs.variance A vector containing the spatially-varying variance
#' corresponding to each data location.
#' @param obs.nuggets A vector containing the spatially-varying nuggets
#' corresponding to each data location.
#'
#' @return This function returns another function for use in \code{optim}.
#'
#' @examples
#' \dontrun{
#' make_global_loglik4_kappa( data, Xmat, cov.model, Scalemat, Distmat, obs.variance, obs.nuggets )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial

make_global_loglik4_kappa <- function( data, Xmat, cov.model, Scalemat, Distmat, obs.variance, obs.nuggets ){

  fixed <- c(FALSE)
  params <- fixed
  function(p){
    params[!fixed] <- p

    # Parameters
    kapp <- params[1]

    Unscl.corr <- cov.spatial( Distmat, cov.model = cov.model,
                               cov.pars = c(1,1), kappa = kapp )
    Corr <- Scalemat*Unscl.corr

    N <- dim(data)[1]
    p <- dim(data)[2]

    Edecomp <- eigen(diag( sqrt(obs.variance) ) %*% Corr %*% diag( sqrt(obs.variance) ) + diag(obs.nuggets))
    Dmat.temp <- diag(Edecomp$values)
    Vmat.temp <- Edecomp$vectors

    Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

    XCX <- t(Xmat) %*% Cov.inv %*% Xmat
    XCX.inv <- chol2inv( chol(XCX) )

    Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

    # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
    loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(data)%*%Ptemp%*%data)) )

    if(abs(loglikelihood) == Inf){ loglikelihood <- 100000 }
    return(loglikelihood)

  }
}

#======================================================================================

#======================================================================================
# Function to calculate the locally anisotropic model
#======================================================================================
# Using a subset of the data, calculate MLEs of covariance and mean
# parameters. Options for estimating models with and without kappa.
#======================================================================================

# Likelihood for Gaussian data with anisotropic covariance,
# without kappa:

#ROxygen comments ----
#' Constructor functions for local parameter estimation.
#'
#' This function generates another function to be used within \code{optim} to
#' obtain maximum likelihood estimates of
#' local variance parameters lam1, lam2, eta, tausq, and sigmasq,
#' assuming the smoothness is fixed, using a Gaussian likelihood with
#' an anisotropic covariance structure.
#'
#' @param locations A matrix of locations.
#' @param cov.model String; the covariance model.
#' @param data A vector or matrix of data to use in the likelihood
#' calculation.
#' @param Xmat The design matrix for the mean model.
#'
#' @return This function returns another function for use in \code{optim}.
#'
#' @examples
#' \dontrun{
#' make_aniso_loglik( locations, cov.model, data, Xmat )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial
#' @importFrom StatMatch mahalanobis.dist

make_aniso_loglik <- function( locations, cov.model, data, Xmat ){

  fixed <- c(FALSE, FALSE, FALSE, FALSE, FALSE)
  params <- fixed
  function(p){
    params[!fixed] <- p

    # Kernel parameters
    lam1 <- params[1]
    lam2 <- params[2]
    eta <- params[3]

    # Nugget
    tausq <- params[4]

    # Process variance
    sigmasq <- params[5]

    # Smoothness
    KAPPA <- NULL

    #================================
    N <- dim(locations)[1]
    p <- dim(data)[2]

    Pmat <- matrix(c(cos(eta),-sin(eta),sin(eta),cos(eta)),nrow=2,byrow=T)
    Dmat <- diag(c(lam1,lam2))

    Sigma <- Pmat %*% Dmat %*% t(Pmat)

    distances <- mahalanobis.dist( data.x = locations, vc = Sigma )
    NS.cov <- sigmasq*cov.spatial(distances, cov.model = cov.model,
                                  cov.pars = c(1,1), kappa = KAPPA)

    Edecomp <- eigen(NS.cov + diag(rep(tausq,N)))
    Dmat.temp <- diag(Edecomp$values)
    Vmat.temp <- Edecomp$vectors

    Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

    XCX <- t(Xmat) %*% Cov.inv %*% Xmat
    XCX.inv <- chol2inv( chol(XCX) )

    Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

    # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
    loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(data)%*%Ptemp%*%data)) )

    if(abs(loglikelihood) == Inf){ loglikelihood <- 1000000 }

    return(loglikelihood)

  }
}


# with kappa:
#ROxygen comments ----
#' Constructor functions for local parameter estimation.
#'
#' This function generates another function to be used within \code{optim} to
#' obtain maximum likelihood estimates of
#' local variance parameters lam1, lam2, eta, tausq, sigmasq, and
#' nu (smoothness) using a Gaussian likelihood with
#' an anisotropic covariance structure.
#'
#' @param locations A matrix of locations.
#' @param cov.model String; the covariance model.
#' @param data A vector or matrix of data to use in the likelihood
#' calculation.
#' @param Xmat The design matrix for the mean model.
#'
#' @return This function returns another function for use in \code{optim}.
#'
#' @examples
#' \dontrun{
#' make_aniso_loglik_kappa( locations, cov.model, data, Xmat )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial
#' @importFrom StatMatch mahalanobis.dist

make_aniso_loglik_kappa <- function( locations, cov.model, data, Xmat ){

  fixed <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
  params <- fixed
  function(p){
    params[!fixed] <- p

    # Kernel parameters
    lam1 <- params[1]
    lam2 <- params[2]
    eta <- params[3]

    # Nugget
    tausq <- params[4]

    # Process variance
    sigmasq <- params[5]

    # Smoothness
    KAPPA <- params[6]

    #================================
    N <- dim(locations)[1]
    p <- dim(data)[2]

    Pmat <- matrix(c(cos(eta),-sin(eta),sin(eta),cos(eta)),nrow=2,byrow=T)
    Dmat <- diag(c(lam1,lam2))

    Sigma <- Pmat %*% Dmat %*% t(Pmat)

    distances <- mahalanobis.dist( data.x = locations, vc = Sigma )
    NS.cov <- sigmasq*cov.spatial(distances, cov.model = cov.model,
                                  cov.pars = c(1,1), kappa = KAPPA)

    Edecomp <- eigen(NS.cov + diag(rep(tausq,N)))
    Dmat.temp <- diag(Edecomp$values)
    Vmat.temp <- Edecomp$vectors

    Cov.inv <- Vmat.temp %*% diag( 1/diag(Dmat.temp) ) %*% t(Vmat.temp)

    XCX <- t(Xmat) %*% Cov.inv %*% Xmat
    XCX.inv <- chol2inv( chol(XCX) )

    Ptemp <- Cov.inv - Cov.inv %*% Xmat %*% XCX.inv %*% t(Xmat) %*% Cov.inv

    # Default is to minimize -- want the maximum. Calculate the neg loglikelihood.
    loglikelihood <- 0.5*( p*sum(log(diag(Dmat.temp))) + p*log(det(XCX)) + sum(diag(t(data)%*%Ptemp%*%data)) )

    if(abs(loglikelihood) == Inf){ loglikelihood <- 1000000 }

    return(loglikelihood)

  }
}

