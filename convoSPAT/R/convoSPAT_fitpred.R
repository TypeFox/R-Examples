#======================================================================================
# Local likelihood estimation for covariance functions with spatially-varying
# parameters: the convoSPAT() package for R
# Mark D. Risser / The Ohio State University / 2014-2015
#======================================================================================

#======================================================================================
# Fit/Predict Functions
#======================================================================================

#======================================================================================
# Fit the nonstationary model
#======================================================================================
# The NSconvo.fit() function estimates the parameters of the nonstationary
# convolution-based spatial model. Required inputs are the observed data and
# locations (a geoR object with $coords and $data).
# Optional inputs include mixture component locations (if not provided, the number of mixture component
# locations are required), the fit radius, the covariance model (exponential is
# the default), and whether or not the nugget and process variance
# will be spatially-varying.
# The output of the model is the mixture component locations, mixture component kernels, estimates of
# mu (mean), tausq (nugget variance), sigmasq (process variance), local MLEs,
# the covariance model, and the MLE covariance matrix.
#======================================================================================
#ROxygen comments ----
#' Fit the nonstationary spatial model
#'
#' \code{NSconvo_fit} estimates the parameters of the nonstationary
#' convolution-based spatial model. Required inputs are the observed data and
#' locations (a geoR object with $coords and $data).
#' Optional inputs include mixture component locations (if not provided,
#' the number of mixture component locations are required), the fit radius,
#' the covariance model (exponential is the default), and whether or not the
#' nugget and process variance will be spatially-varying.
#'
#' @param geodata A list containing elements \code{coords} and \code{data} as
#' described next. Typically an object of the class "\code{geodata}", although
#' a geodata object only allows \code{data} to be a vector (no replicates).
#' If not provided, the arguments \code{coords} and \code{data} must be
#' provided instead.
#' @param sp.SPDF A "\code{SpatialPointsDataFrame}" object, which contains the
#' spatial coordinates and additional attribute variables corresponding to the
#' spatoal coordinates
#' @param coords An N x 2 matrix where each row has the two-dimensional
#' coordinates of the N data locations. By default, it takes the \code{coords}
#' component of the argument \code{geodata}, if provided.
#' @param data A vector or matrix with N rows, containing the data values.
#' Inputting a vector corresponds to a single replicate of data, while
#' inputting a matrix corresponds to replicates. In the case of replicates,
#' the model assumes the replicates are independent and identically
#' distributed.
#' @param cov.model A string specifying the model for the correlation
#' function; following \code{geoR}, defaults to \code{"exponential"}.
#' Options available in this package are: "\code{exponential}",
#' \code{"cauchy"}, \code{"matern"}, \code{"circular"}, \code{"cubic"},
#' \code{"gaussian"}, \code{"spherical"}, and \code{"wave"}. For further
#' details, see documentation for \code{\link[geoR]{cov.spatial}}.
#' @param mean.model An object of class \code{\link[stats]{formula}},
#' specifying the mean model to be used. Defaults to an intercept only.
#' @param mc.locations Optional; matrix of mixture component locations.
#' @param N.mc Optional; if \code{mc.locations} is not specified, the
#' function will create a rectangular grid of size \code{N.mc} over the
#' spatial domain.
#' @param lambda.w Scalar; tuning parameter for the weight function.
#' Defaults to be the square of one-half of the minimum distance between
#' mixture component locations.
#' @param mc.kernels Optional specification of mixture component kernel
#' matrices (based on expert opinion, etc.).
#' @param fit.radius Scalar; specifies the fit radius or neighborhood size
#' for the local likelihood estimation.
#' @param ns.nugget Logical; indicates if the nugget variance (tausq) should
#' be spatially-varying (\code{TRUE}) or constant (\code{FALSE}).
#' @param ns.variance Logical; indicates if the process variance (sigmasq)
#' should be spatially-varying (\code{TRUE}) or constant (\code{FALSE}).
#'
#' @param local.pars.LB,local.pars.UB Optional vectors of lower and upper
#' bounds, respectively, used by the \code{"L-BFGS-B"} method option in the
#' \code{\link[stats]{optim}} function for the local parameter estimation.
#' Each vector must be of length five,
#' containing values for lam1, lam2, tausq, sigmasq, and nu. Default for
#' \code{local.pars.LB} is \code{rep(1e-05,5)}; default for
#' \code{local.pars.UB} is \code{c(max.distance/2, max.distance/2, 4*resid.var, 4*resid.var, 100)},
#' where \code{max.distance} is the maximum interpoint distance of the
#' observed data and \code{resid.var} is the residual variance from using
#' \code{\link[stats]{lm}} with \code{mean.model}.
#'
#' @param global.pars.LB,global.pars.UB Optional vectors of lower and upper
#' bounds, respectively, used by the \code{"L-BFGS-B"} method option in the
#' \code{\link[stats]{optim}} function for the global parameter estimation.
#' Each vector must be of length three,
#' containing values for tausq, sigmasq, and nu. Default for
#' \code{global.pars.LB} is \code{rep(1e-05,3)}; default for
#' \code{global.pars.UB} is \code{c(4*resid.var, 4*resid.var, 100)},
#' where \code{resid.var} is the residual variance from using
#' \code{\link[stats]{lm}} with \code{mean.model}.
#'
#' @param local.ini.pars Optional vector of initial values used by the
#' \code{"L-BFGS-B"} method option in the \code{\link[stats]{optim}}
#' function for the local parameter estimation. The vector must be of length
#' five, containing values for lam1, lam2, tausq, sigmasq, and nu. Defaults
#' to \code{c(max.distance/10, max.distance/10, 0.1*resid.var, 0.9*resid.var, 1)},
#' where \code{max.distance} is the maximum interpoint distance of the
#' observed data and \code{resid.var} is the residual variance from using
#' \code{\link[stats]{lm}} with \code{mean.model}.
#'
#' @param global.ini.pars Optional vector of initial values used by the
#' \code{"L-BFGS-B"} method option in the \code{\link[stats]{optim}}
#' function for the global parameter estimation. The vector must be of length
#' three, containing values for tausq, sigmasq, and nu. Defaults to
#' \code{c(0.1*resid.var, 0.9*resid.var, 1)}, where \code{resid.var} is the
#' residual variance from using \code{\link[stats]{lm}} with \code{mean.model}.
#'
#' @return A "NSconvo" object, with the following components:
#' \item{mc.locations}{Mixture component locations used for the simulated
#' data.}
#' \item{mc.kernels}{Mixture component kernel matrices used for the simulated
#' data.}
#' \item{MLEs.save}{Table of local maximum likelihood estimates for each
#' mixture component location.}
#' \item{kernel.ellipses}{\code{N.obs} x 2 x 2 array, containing the kernel
#' matrices corresponding to each of the simulated values.}
#' \item{data}{Observed data values.}
#' \item{beta.GLS}{Vector of generalized least squares estimates of beta,
#' the mean coefficients.}
#' \item{beta.cov}{Covariance matrix of the generalized least squares
#' estimate of beta.}
#' \item{Mean.coefs}{"Regression table" for the mean coefficient estimates,
#' listing the estimate, standard error, and t-value.}
#' \item{tausq.est}{Estimate of tausq (nugget variance), either scalar (when
#' \code{ns.nugget = "FALSE"}) or a vector of length N (when
#' \code{ns.nugget = "TRUE"}), which contains the estimated nugget variance
#' for each observation location.}
#' \item{sigmasq.est}{Estimate of sigmasq (process variance), either scalar
#' (when \code{ns.variance = "FALSE"}) or a vector of length N (when
#' \code{ns.variance = "TRUE"}), which contains the estimated process
#' variance for each observation location.}
#' \item{kappa.MLE}{Scalar maximum likelihood estimate for kappa (when
#' applicable).}
#' \item{Cov.mat}{Estimated covariance matrix (\code{N.obs} x \code{N.obs})
#' using all relevant parameter estimates.}
#' \item{Cov.mat.inv}{Inverse of \code{Cov.mat}, the estimated covariance
#' matrix (\code{N.obs} x \code{N.obs}).}
#' \item{cov.model}{String; the correlation model used for estimation.}
#' \item{ns.nugget}{Logical, indicating if the nugget variance was estimated
#' as spatially-varing (\code{TRUE}) or constant (\code{FALSE}).}
#' \item{ns.variance}{Logical, indicating if the process variance was
#' estimated as spatially-varying (\code{TRUE}) or constant (\code{FALSE}).}
#' \item{coords}{N x 2 matrix of observation locations.}
#' \item{global.loglik}{Scalar value of the maximized likelihood from the
#' global optimization (if available).}
#' \item{Xmat}{Design matrix, obtained from using \code{\link[stats]{lm}}
#' with \code{mean.model}.}
#' \item{lambda.w}{Tuning parameter for the weight function.}
#'
#' @examples
#' # Using white noise data
#' fit.model <- NSconvo_fit( coords = cbind( runif(100), runif(100)),
#' data = rnorm(100), fit.radius = 0.4, N.mc = 4 )
#'
#' @export
#' @importFrom geoR cov.spatial
#' @importFrom stats lm
#' @importFrom stats optim

NSconvo_fit <- function( geodata = NULL, sp.SPDF = NULL,
                         coords = geodata$coords, data = geodata$data,
                         cov.model = "exponential", mean.model = data ~ 1,
                         mc.locations = NULL, N.mc = NULL, lambda.w = NULL,
                         mc.kernels = NULL, fit.radius = NULL,
                         ns.nugget = FALSE, ns.variance = FALSE,
                         local.pars.LB = NULL, local.pars.UB = NULL,
                         global.pars.LB = NULL, global.pars.UB = NULL,
                         local.ini.pars = NULL, global.ini.pars = NULL ){

  if( is.null(fit.radius) ){cat("\nPlease specify a fitting radius.\n")}
  #===========================================================================
  # Formatting for coordinates/data
  #===========================================================================
  if( is.null(geodata) == FALSE ){
    if( class(geodata) != "geodata" ){
      cat("\nPlease use a geodata object for the 'geodata = ' input.\n")
    }
    coords <- geodata$coords
    data <- geodata$data
  }
  if( is.null(sp.SPDF) == FALSE ){
    if( class(sp.SPDF) != "SpatialPointsDataFrame" ){
      cat("\nPlease use a SpatialPointsDataFrame object for the 'sp.SPDF = ' input.\n")
    }
    geodata <- geoR::as.geodata( sp.SPDF )
    coords <- geodata$coords
    data <- geodata$data
  }

  coords <- as.matrix(coords)
  N <- dim(coords)[1]
  data <- as.matrix(data, nrow=N)
  p <- dim(data)[2]

  # Make sure cov.model is one of the permissible options
  if( cov.model != "cauchy" & cov.model != "matern" & cov.model != "circular" &
        cov.model != "cubic" & cov.model != "gaussian" & cov.model != "exponential" &
        cov.model != "spherical" & cov.model != "wave" ){
    cat("\nPlease specify a valid covariance model (cauchy, matern,\ncircular, cubic, gaussian,
          exponential, spherical, or wave).\n")
  }

  #===========================================================================
  # Calculate the mixture component locations if not user-specified
  #===========================================================================
  if( is.null(mc.locations) == TRUE ){ # Calculate the mixture component locations

    if( is.null(N.mc) == TRUE ){
      cat("Please enter the desired number of mixture component locations. \n")
    }

    lon_min <- min(coords[,1])
    lon_max <- max(coords[,1])
    lat_min <- min(coords[,2])
    lat_max <- max(coords[,2])

    #=======================================
    # mixture component knot locations
    #=======================================
    mc_x <- seq(from = lon_min + 0.5*(lon_max - lon_min)/floor(sqrt(N.mc)),
                   to = lon_max - 0.5*(lon_max - lon_min)/floor(sqrt(N.mc)),
                   length = floor(sqrt(N.mc)) )
    mc_y <- seq(from = lat_min + 0.5*(lat_max - lat_min)/floor(sqrt(N.mc)),
                   to = lat_max - 0.5*(lat_max - lat_min)/floor(sqrt(N.mc)),
                   length = floor(sqrt(N.mc)) )
    mc.locations <- expand.grid( mc_x, mc_y )
    mc.locations <- matrix(c(mc.locations[,1], mc.locations[,2]), ncol=2, byrow=F)
  }

  K <- dim(mc.locations)[1]

  #===========================================================================
  # Check the mixture component locations
  #===========================================================================
  check.mc.locs <- mc_N( coords, mc.locations, fit.radius )
  cat(paste("The", K, "local models will be fit with local sample sizes\nranging between", min(check.mc.locs),
            "and", max(check.mc.locs), ".\n" ))
  if( min(check.mc.locs) < 5 ){cat("\nWARNING: at least one of the mc locations has too few data points.\n")}

  #===========================================================================
  # Set the tuning parameter, if not specified
  #===========================================================================
  if( is.null(lambda.w) == TRUE ){
    lambda.w <- ( 0.5*sqrt(sum((mc.locations[1,] - mc.locations[2,])^2)) )^2
  }

  #===========================================================================
  # Calculate the design matrix
  #===========================================================================
  OLS.model <- lm( mean.model, x=TRUE )

  Xmat <- matrix( unname( OLS.model$x ), nrow=N )

  #===========================================================================
  # Specify lower, upper, and initial parameter values for optim()
  #===========================================================================
  lon_min <- min(coords[,1])
  lon_max <- max(coords[,1])
  lat_min <- min(coords[,2])
  lat_max <- max(coords[,2])

  max.distance <- sqrt( sum((c(lon_min,lat_min) - c(lon_max,lat_max))^2))

  if( p > 1 ){
    ols.sigma <- NULL
    for(i in 1:length(names(summary(OLS.model)))){
      ols.sigma <- c( ols.sigma, summary(OLS.model)[[i]]$sigma )
    }
    resid.var <- (max(ols.sigma))^2
  }
  if( p == 1 ){
    resid.var <- summary(OLS.model)$sigma^2
  }

  #=================================
  # Lower limits for optim()
  #=================================
  if( is.null(local.pars.LB) == TRUE ){
    lam1.LB <- 1e-05
    lam2.LB <- 1e-05
    tausq.local.LB <- 1e-05
    sigmasq.local.LB <- 1e-05
    kappa.local.LB <- 1e-05
  }
  if( is.null(local.pars.LB) == FALSE ){
    lam1.LB <- local.pars.LB[1]
    lam2.LB <- local.pars.LB[2]
    tausq.local.LB <- local.pars.LB[3]
    sigmasq.local.LB <- local.pars.LB[4]
    kappa.local.LB <- local.pars.LB[5]
  }
  if( is.null(global.pars.LB) == TRUE ){
    tausq.global.LB <- 1e-05
    sigmasq.global.LB <- 1e-05
    kappa.global.LB <- 1e-05
  }
  if( is.null(global.pars.LB) == FALSE ){
    tausq.global.LB <- global.pars.LB[1]
    sigmasq.global.LB <- global.pars.LB[2]
    kappa.global.LB <- global.pars.LB[3]
  }

  #=================================
  # Upper limits for optim()
  #=================================
  if( is.null(local.pars.UB) == TRUE ){
    lam1.UB <- max.distance/4
    lam2.UB <- max.distance/4
    tausq.local.UB <- 4*resid.var
    sigmasq.local.UB <- 4*resid.var
    kappa.local.UB <- 30
  }
  if( is.null(local.pars.UB) == FALSE ){
    lam1.UB <- local.pars.UB[1]
    lam2.UB <- local.pars.UB[2]
    tausq.local.UB <- local.pars.UB[3]
    sigmasq.local.UB <- local.pars.UB[4]
    kappa.local.UB <- local.pars.UB[5]
  }
  if( is.null(global.pars.UB) == TRUE ){
    tausq.global.UB <- 4*resid.var
    sigmasq.global.UB <- 4*resid.var
    kappa.global.UB <- 30
  }
  if( is.null(global.pars.UB) == FALSE ){
    tausq.global.UB <- global.pars.UB[1]
    sigmasq.global.UB <- global.pars.UB[2]
    kappa.global.UB <- global.pars.UB[3]
  }

  #=================================
  # Initial values for optim()
  #=================================

  # Local estimation
  if( is.null(local.ini.pars) == TRUE ){
    lam1.init <- max.distance/10
    lam2.init <- max.distance/10
    tausq.local.init <- 0.1*resid.var
    sigmasq.local.init <- 0.9*resid.var
    kappa.local.init <- 1
  }
  if( is.null(local.ini.pars) == FALSE ){
    lam1.init <- local.ini.pars[1]
    lam2.init <- local.ini.pars[2]
    tausq.local.init <- local.ini.pars[3]
    sigmasq.local.init <- local.ini.pars[4]
    kappa.local.init <- local.ini.pars[5]
  }

  # Global estimation
  if( is.null(global.ini.pars) == TRUE ){
    tausq.global.init <- 0.1*resid.var
    sigmasq.global.init <- 0.9*resid.var
    kappa.global.init <- 1
  }
  if( is.null(global.ini.pars) == FALSE ){
    tausq.global.init <- global.ini.pars[1]
    sigmasq.global.init <- global.ini.pars[2]
    kappa.global.init <- global.ini.pars[3]
  }

  #===========================================================================
  # Calculate the mixture component kernels if not user-specified
  #===========================================================================
  if( is.null(mc.kernels) == TRUE ){

    # Storage for the mixture component kernels
    mc.kernels <- array(NA, dim=c(2,2,K))
    MLEs.save <- matrix(NA,K,7)

    # Estimate the kernel function for each mixture component location,
    # completely specified by the kernel covariance matrix
    for( k in 1:K ){

      temp.locs <- coords[ abs(coords[,1]-mc.locations[k,1]) <= fit.radius
                           & (abs(coords[,2] - mc.locations[k,2]) <= fit.radius), ]
      temp.dat <- data[(abs(coords[,1] - mc.locations[k,1]) <= fit.radius)
                       & (abs(coords[,2] - mc.locations[k,2]) <= fit.radius), ]
      X.tem <- as.matrix(Xmat[ abs(coords[,1]-mc.locations[k,1]) <= fit.radius
                               & (abs(coords[,2] - mc.locations[k,2]) <= fit.radius), ])

      # Isolate the data/locations to be used for calculating the local kernel
      distances <- rep(NA,dim(temp.locs)[1])

      for(i in 1:dim(temp.locs)[1]){
        distances[i] <- sqrt(sum((temp.locs[i,] - mc.locations[k,])^2))
      }

      temp.locations <- temp.locs[distances <= fit.radius,]
      Xtemp <- X.tem[distances <= fit.radius,]
      n.fit <- dim(temp.locations)[1]
      temp.dat <- as.matrix( temp.dat, nrow=n.fit)
      temp.data <- temp.dat[distances <= fit.radius,]
      temp.data <- as.matrix(temp.data, nrow=n.fit)

      cat("Calculating the mixture component parameter set for location ", k," using ",
                  n.fit," observations. \n", sep="")

      # Covariance models with the kappa parameter
      if( cov.model == "matern" || cov.model == "cauchy" ){

        # Calculate a locally anisotropic covariance
        # Parameter order is lam1, lam2, eta, tausq, sigmasq, kappa
        anis.model.kappa <- make_aniso_loglik_kappa( locations = temp.locations,
                                                     cov.model = cov.model,
                                                     data = temp.data,
                                                     Xmat = Xtemp )

        MLEs <- optim( c(lam1.init, lam2.init, pi/4, tausq.local.init,
                         sigmasq.local.init, kappa.local.init ),
                       anis.model.kappa,
                       method = "L-BFGS-B",
                       lower=c( lam1.LB, lam2.LB, 0, tausq.local.LB, sigmasq.local.LB ),
                       upper=c( lam1.UB, lam2.UB, pi/2,
                                tausq.local.UB, sigmasq.local.UB, kappa.local.UB ) )
        if( MLEs$convergence != 0 ){
          cat( "There was an error with optim(): \n", MLEs$message, "\n" )
        }
        MLE.pars <- MLEs$par
      }
      # Covariance models without the kappa parameter
      if( cov.model != "matern" & cov.model != "cauchy" ){

        # Calculate a locally anisotropic covariance
        # Parameter order is lam1, lam2, eta, tausq, sigmasq
        anis.model <- make_aniso_loglik( locations = temp.locations,
                                         cov.model = cov.model,
                                         data = temp.data,
                                         Xmat = Xtemp )

        MLEs <- optim( c( lam1.init, lam2.init, pi/4, tausq.local.init,
                          sigmasq.local.init ),
                       anis.model,
                       method = "L-BFGS-B",
                       lower=c( lam1.LB, lam2.LB, 0, tausq.local.LB, sigmasq.local.LB ),
                       upper=c( lam1.UB, lam2.UB, pi/2,
                                tausq.local.UB, sigmasq.local.UB ) )
        if( MLEs$convergence != 0 ){
          cat( "There was an error with optim(): \n", MLEs$message, "\n" )
        }

        MLE.pars <- c(MLEs$par, NA)
      }

      # Save the kernel matrix
      mc.kernels[,,k] <- kernel_cov( MLE.pars[1:3] )

      # Save all MLEs
      # Parameter order: n, lam1, lam2, eta, tausq, sigmasq, mu, kappa
      MLEs.save[k,] <- c( n.fit, MLE.pars )
    }

    # Put the MLEs into a data frame
    MLEs.save <- data.frame( MLEs.save )
    names(MLEs.save) <- c("n","lam1", "lam2", "eta", "tausq", "sigmasq", "kappa" )
  }

  #===========================================================================
  # Calculate the weights for each observation location
  #===========================================================================
  weights <- matrix(NA, N, K)
  for(n in 1:N){for(k in 1:K){
    weights[n,k] <- exp(-sum((coords[n,] - mc.locations[k,])^2)/(2*lambda.w))
  }
                # Normalize the weights
                weights[n,] <- weights[n,]/sum(weights[n,])
  }

  #===========================================================================
  # Calculate the kernel ellipses
  #===========================================================================
  kernel.ellipses <- array(0, dim=c(2,2,N))
  for(n in 1:N){
    for(k in 1:K){
      kernel.ellipses[,,n] <- kernel.ellipses[,,n] + weights[n,k]*mc.kernels[,,k]
    }
  }

  #===============
  # If specified: calculate the spatially-varying nugget and variance
  if( ns.nugget == TRUE   ){
    mc.nuggets <- as.numeric(MLEs.save$tausq)
    obs.nuggets <- rep(0,N)
    for(n in 1:N){
      for(k in 1:K){
        obs.nuggets[n] <- obs.nuggets[n] + weights[n,k]*mc.nuggets[k]
      }
    }
  }

  if( ns.variance == TRUE ){
    mc.variance <- as.numeric(MLEs.save$sigmasq)
    obs.variance <- rep(0,N)
    for(n in 1:N){
      for(k in 1:K){
        obs.variance[n] <- obs.variance[n] + weights[n,k]*mc.variance[k]
      }
    }
  }

  #===========================================================================
  # Global parameter estimation
  #===========================================================================
  # First, for models without kappa.
  if( cov.model != "matern" & cov.model != "cauchy" ){

    KAPPA <- NULL
    #====================================
    # Calculate the correlation matrix

    Scale.mat <- matrix(rep(NA, N^2), nrow=N)
    Dist.mat <- matrix(rep(NA, N^2), nrow=N)
    # Calculate the elements of the observed correlations matrix.
    for(i in 1:N){

      # Diagonal elements
      Kerneli <- kernel.ellipses[,,i]
      det_i <- Kerneli[1,1]*Kerneli[2,2] - Kerneli[1,2]*Kerneli[2,1]

      Scale.mat[i,i] <- 1
      Dist.mat[i,i] <- 0

      # Ui <- chol(Kerneli)

      if(i < N){
        for(j in (i+1):N){ # Off-diagonal elements

          Kernelj <- kernel.ellipses[,,j]
          det_j <- Kernelj[1,1]*Kernelj[2,2] - Kernelj[1,2]*Kernelj[2,1]

          avg_ij <- 0.5 * (Kerneli + Kernelj)
          # Uij <- chol(avg_ij)
          det_ij <- avg_ij[1,1]*avg_ij[2,2] - avg_ij[1,2]*avg_ij[2,1]
          #vec_ij <- backsolve(Uij, (coords[i,]-coords[j,]), transpose = TRUE)

          Scale.mat[i,j] <- sqrt( sqrt(det_i*det_j) / det_ij )
          # Dist.mat[i,j] <- sqrt(sum(vec_ij^2))
          Dist.mat[i,j] <- sqrt( t(coords[i,]-coords[j,]) %*% solve(avg_ij) %*% (coords[i,]-coords[j,]) )

          Scale.mat[j,i] <- Scale.mat[i,j]
          Dist.mat[j,i] <- Dist.mat[i,j]

        }
      }
    }
    Unscl.corr <- cov.spatial( Dist.mat, cov.model = cov.model,
                               cov.pars = c(1,1), kappa = KAPPA )
    NS.corr <- Scale.mat*Unscl.corr

    #====================================
    # Global parameter estimation

    if( ns.nugget == FALSE & ns.variance == FALSE ){
      cat("Calculating the variance parameter MLEs. \n")

      Corr.decomp <- eigen(NS.corr)

      Dmat <- diag(Corr.decomp$values)
      Vmat <- Corr.decomp$vectors

      overall.lik1 <- make_global_loglik1( data = data, Xmat = Xmat,
                                           Dmat = Dmat, Vmat = Vmat )

      overall.MLEs <- optim( c( tausq.global.init, sigmasq.global.init ),
                             overall.lik1,
                             method = "L-BFGS-B",
                             lower=c( tausq.global.LB, sigmasq.global.LB ),
                             upper=c( tausq.global.UB, sigmasq.global.UB ) )
      if( overall.MLEs$convergence != 0 ){
        cat( "There was an error with optim(): \n", overall.MLEs$message, "\n" )
      }
      tausq.MLE <- overall.MLEs$par[1]
      sigmasq.MLE <- overall.MLEs$par[2]

      global.lik <- overall.MLEs$value

      ObsNuggMat <- diag(rep(tausq.MLE,N))
      ObsCov <- sigmasq.MLE * NS.corr

      obs.variance <- rep(sigmasq.MLE, N)
    }

    if( ns.nugget == TRUE & ns.variance == FALSE ){
      cat("Calculating the variance parameter MLEs. \n")

      overall.lik2 <- make_global_loglik2( data = data,
                                           Xmat = Xmat,
                                           Corr = NS.corr,
                                           obs.nuggets = obs.nuggets )

      overall.MLEs <- optim( sigmasq.global.init, overall.lik2, method = "L-BFGS-B",
                             lower=c( sigmasq.global.LB ),
                             upper=c( sigmasq.global.UB ) )
      if( overall.MLEs$convergence != 0 ){
        cat( "There was an error with optim(): \n", overall.MLEs$message, "\n" )
      }
      sigmasq.MLE <- overall.MLEs$par[1]
      global.lik <- overall.MLEs$value

      ObsNuggMat <- diag(obs.nuggets)
      ObsCov <- sigmasq.MLE * NS.corr

      obs.variance <- rep(sigmasq.MLE, N)

    }

    if( ns.nugget == FALSE & ns.variance == TRUE ){
      cat("Calculating the variance parameter MLEs. \n")

      overall.lik3 <- make_global_loglik3( data = data,
                                           Xmat = Xmat,
                                           Corr = NS.corr,
                                           obs.variance = obs.variance )

      overall.MLEs <- optim( tausq.global.init, overall.lik3, method = "L-BFGS-B",
                             lower=c( tausq.global.LB ),
                             upper=c( tausq.global.UB ) )

      if( overall.MLEs$convergence != 0 ){
        cat( "There was an error with optim(): \n", overall.MLEs$message, "\n" )
      }
      tausq.MLE <- overall.MLEs$par[1]
      global.lik <- overall.MLEs$value

      ObsNuggMat <- diag(rep(tausq.MLE,N))
      ObsCov <- diag(sqrt(obs.variance)) %*% NS.corr %*% diag(sqrt(obs.variance))

    }

    if( ns.nugget == TRUE & ns.variance == TRUE ){

      Cov <- diag( sqrt(obs.variance) ) %*% NS.corr %*% diag( sqrt(obs.variance) )

      global.lik <- NA
      ObsNuggMat <- diag(obs.nuggets)
      ObsCov <- Cov

    }
    kappa.MLE <- NA
  }

  # Next, for models with kappa.
  if( cov.model == "matern" || cov.model == "cauchy" ){

    #====================================
    # Calculate the correlation matrix

    Scale.mat <- matrix(rep(NA, N^2), nrow=N)
    Dist.mat <- matrix(rep(NA, N^2), nrow=N)
    # Calculate the elements of the observed correlations matrix.
    for(i in 1:N){

      # Diagonal elements
      Kerneli <- kernel.ellipses[,,i]
      det_i <- Kerneli[1,1]*Kerneli[2,2] - Kerneli[1,2]*Kerneli[2,1]

      Scale.mat[i,i] <- 1
      Dist.mat[i,i] <- 0

      #Ui <- chol(Kerneli)

      if(i < N){
        for(j in (i+1):N){ # Off-diagonal elements

          Kernelj <- kernel.ellipses[,,j]
          det_j <- Kernelj[1,1]*Kernelj[2,2] - Kernelj[1,2]*Kernelj[2,1]

          avg_ij <- 0.5 * (Kerneli + Kernelj)
          # Uij <- chol(avg_ij)
          det_ij <- avg_ij[1,1]*avg_ij[2,2] - avg_ij[1,2]*avg_ij[2,1]
          #vec_ij <- backsolve(Uij, (coords[i,]-coords[j,]), transpose = TRUE)

          Scale.mat[i,j] <- sqrt( sqrt(det_i*det_j) / det_ij )
          # Dist.mat[i,j] <- sqrt(sum(vec_ij^2))
          Dist.mat[i,j] <- sqrt( t(coords[i,]-coords[j,]) %*% solve(avg_ij) %*% (coords[i,]-coords[j,]) )

          Scale.mat[j,i] <- Scale.mat[i,j]
          Dist.mat[j,i] <- Dist.mat[i,j]

        }
      }
    }

    #====================================
    # Global parameter estimation

    if( ns.nugget == FALSE & ns.variance == FALSE ){
      cat("Calculating the variance and smoothness parameter MLEs. \n")

      overall.lik1.kappa <- make_global_loglik1_kappa( data = data, Xmat = Xmat,
                                                       cov.model = cov.model,
                                                       Scalemat = Scale.mat,
                                                       Distmat = Dist.mat )

      overall.MLEs <- optim(c( tausq.global.init, sigmasq.global.init, kappa.global.init ),
                            overall.lik1.kappa,
                            method = "L-BFGS-B",
                            lower=c( tausq.global.LB, sigmasq.global.LB, kappa.global.LB ),
                            upper=c( tausq.global.UB, sigmasq.global.UB, kappa.global.UB ) )
      if( overall.MLEs$convergence != 0 ){
        cat( "There was an error with optim(): \n", overall.MLEs$message, "\n" )
      }
      tausq.MLE <- overall.MLEs$par[1]
      sigmasq.MLE <- overall.MLEs$par[2]
      kappa.MLE <- overall.MLEs$par[3]

      global.lik <- overall.MLEs$value
      ObsNuggMat <- diag(rep(tausq.MLE,N))
      ObsCov <- sigmasq.MLE * Scale.mat * cov.spatial( Dist.mat, cov.model = cov.model,
                                                       cov.pars = c(1,1), kappa = kappa.MLE )
      obs.variance <- rep(sigmasq.MLE, N)
    }

    if( ns.nugget == TRUE & ns.variance == FALSE ){
      cat("Calculating the variance and smoothness parameter MLEs. \n")

      overall.lik2.kappa <- make_global_loglik2_kappa( data = data, Xmat = Xmat,
                                                       cov.model = cov.model,
                                                       Scalemat = Scale.mat,
                                                       Distmat = Dist.mat,
                                                       obs.nuggets = obs.nuggets )

      overall.MLEs <- optim(c( sigmasq.global.init, kappa.global.init ),
                            overall.lik2.kappa,
                            method = "L-BFGS-B",
                            lower=c( sigmasq.global.LB, kappa.global.LB ),
                            upper=c( sigmasq.global.UB, kappa.global.UB ) )
      if( overall.MLEs$convergence != 0 ){
        cat( "There was an error with optim(): \n", overall.MLEs$message, "\n" )
      }
      sigmasq.MLE <- overall.MLEs$par[1]
      kappa.MLE <- overall.MLEs$par[2]

      global.lik <- overall.MLEs$value
      ObsNuggMat <- diag(obs.nuggets)
      ObsCov <- sigmasq.MLE * Scale.mat * cov.spatial( Dist.mat, cov.model = cov.model,
                                                       cov.pars = c(1,1), kappa = kappa.MLE )

      obs.variance <- rep(sigmasq.MLE, N)

    }

    if( ns.nugget == FALSE & ns.variance == TRUE ){
      cat("Calculating the variance and smoothness parameter MLEs. \n")

      overall.lik3.kappa <- make_global_loglik3_kappa( data = data, Xmat = Xmat,
                                                       cov.model = cov.model,
                                                       Scalemat = Scale.mat,
                                                       Distmat = Dist.mat,
                                                       obs.variance = obs.variance )

      overall.MLEs <- optim(c( tausq.global.init, kappa.global.init ),
                            overall.lik3.kappa,
                            method = "L-BFGS-B",
                            lower=c( tausq.global.LB, kappa.global.LB ),
                            upper=c( tausq.global.UB, kappa.global.UB ) )
      if( overall.MLEs$convergence != 0 ){
        cat( "There was an error with optim(): \n", overall.MLEs$message, "\n" )
      }
      tausq.MLE <- overall.MLEs$par[1]
      kappa.MLE <- overall.MLEs$par[2]

      global.lik <- overall.MLEs$value
      ObsNuggMat <- diag(rep(tausq.MLE,N))
      ObsCov <- diag(sqrt(obs.variance)) %*% (Scale.mat * cov.spatial( Dist.mat, cov.model = cov.model,
                                                                       cov.pars = c(1,1), kappa = kappa.MLE )) %*% diag(sqrt(obs.variance))

    }

    if( ns.nugget == TRUE & ns.variance == TRUE ){
      cat("Calculating the smoothness parameter MLEs. \n")

      overall.lik4.kappa <- make_global_loglik4_kappa( data = data, Xmat = Xmat,
                                                       cov.model = cov.model,
                                                       Scalemat = Scale.mat,
                                                       Distmat = Dist.mat,
                                                       obs.nuggets = obs.nuggets,
                                                       obs.variance = obs.variance )

      overall.MLEs <- optim( kappa.global.init, overall.lik4.kappa, method = "L-BFGS-B",
                             lower=c( kappa.global.LB ),
                             upper=c( kappa.global.UB ) )
      if( overall.MLEs$convergence != 0 ){
        cat( "There was an error with optim(): \n", overall.MLEs$message, "\n" )
      }

      kappa.MLE <- overall.MLEs$par[1]

      global.lik <- overall.MLEs$value

      ObsNuggMat <- diag(obs.nuggets)
      ObsCov <- diag(sqrt(obs.variance)) %*% (Scale.mat * cov.spatial( Dist.mat, cov.model = cov.model,
                                                                       cov.pars = c(1,1), kappa = kappa.MLE )) %*% diag(sqrt(obs.variance))

    }
  }

  Data.Cov <- ObsNuggMat + ObsCov
  #===========================================================================
  # Calculate the GLS estimate of beta
  #===========================================================================

  Data.Cov.inv <- solve(Data.Cov)
  beta.cov <- chol2inv( chol(t(Xmat)%*%Data.Cov.inv%*%Xmat) )/p
  beta.GLS <- (p*beta.cov %*% t(Xmat) %*% Data.Cov.inv %*% data %*% rep(1,p))/p
  Mean.coefs <- data.frame( Estimate = beta.GLS,
                            Std.Error = sqrt(diag(beta.cov)),
                            t.val = beta.GLS/sqrt(diag(beta.cov)) )

  #===========================================================================
  # Output
  #===========================================================================
  if( ns.nugget == TRUE ){
    tausq.out <- obs.nuggets
  }
  if( ns.nugget == FALSE ){
    tausq.out <- tausq.MLE
  }

  if( ns.variance == TRUE ){
    sigmasq.out <- obs.variance
  }
  if( ns.variance == FALSE ){
    sigmasq.out <- sigmasq.MLE
  }

  output <- list( mc.kernels = mc.kernels,
                  mc.locations = mc.locations,
                  MLEs.save = MLEs.save,
                  kernel.ellipses = kernel.ellipses,
                  data = data,
                  beta.GLS = beta.GLS,
                  beta.cov = beta.cov,
                  Mean.coefs = Mean.coefs,
                  tausq.est = tausq.out,
                  sigmasq.est = sigmasq.out,
                  kappa.MLE = kappa.MLE,
                  Cov.mat = Data.Cov,
                  Cov.mat.inv = Data.Cov.inv,
                  cov.model = cov.model,
                  ns.nugget = ns.nugget,
                  ns.variance = ns.variance,
                  coords = coords,
                  global.loglik = global.lik,
                  Xmat = Xmat,
                  lambda.w = lambda.w )

  class(output) <- "NSconvo"

  return(output)

}

#======================================================================================
# Calculate predictions using the output of NSconvo_fit()
#======================================================================================
# Using the output from NSconvo_fit(), calculate the kriging predictors
# and kriging standard errors for prediction locations of interest.
#======================================================================================
#ROxygen comments ----
#' Obtain predictions at unobserved locations for the nonstationary
#' spatial model.
#'
#' \code{predict.NSconvo} calculates the kriging predictor and corresponding
#' standard errors at unmonitored sites.
#'
#' @param object A "NSconvo" object, from \code{NSconvo_fit}.
#' @param pred.coords Matrix of locations where predictions are required.
#' @param pred.covariates Matrix of covariates for the prediction locations,
#' NOT including an intercept. The number of columns for this matrix must
#' match the design matrix from \code{mean.model} in \code{\link{NSconvo_fit}}.
#' Defaults to an intercept only.
#' @param ... additional arguments affecting the predictions produced.
#'
#' @return A list with the following components:
#' \item{pred.means}{Vector of the kriging predictor, for each location in
#' \code{pred.coords}.}
#' \item{pred.SDs}{Vector of the kriging standard errors, for each location
#' in \code{pred.coords}.}
#'
#' @examples
#' \dontrun{
#' pred.NS <- predict( NSconvo.obj,
#' pred.coords = matrix(c(1,1), ncol=2),
#' pred.covariates = matrix(c(1,1), ncol=2) )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial
#'

predict.NSconvo <- function(object, pred.coords, pred.covariates = NULL,
                            ... )
{
  if( !inherits(object, "NSconvo") ){
    warning("Object is not of type NSconvo.")
  }
  else{
    {
      pred.coords <- as.matrix(pred.coords)
      M <- dim(pred.coords)[1]
      mc.locations <- object$mc.locations
      K <- dim(mc.locations)[1]
      mc.kernels <- object$mc.kernels
      kernel.ellipses <- object$kernel.ellipses
      ns.nugget <- object$ns.nugget
      ns.variance <- object$ns.variance
      beta.MLE <- as.matrix(object$beta.GLS)
      tausq.est <- object$tausq.est
      sigmasq.est <- object$sigmasq.est
      kappa.MLE <- object$kappa.MLE
      mc.MLEs <- object$MLEs.save
      Cov.mat.inv <- object$Cov.mat.inv
      data <- object$data
      N <- length(object$data)
      coords <- object$coords
      cov.model <- object$cov.model
      Xmat <- object$Xmat
      lambda.w <- object$lambda.w
      if (is.null(pred.covariates) == TRUE) {
        Xpred <- rep(1, M)
      }
      if (is.null(pred.covariates) == FALSE) {
        Xpred <- cbind(rep(1, M), pred.covariates)
      }
      data <- matrix(data, nrow = N)
      p <- dim(data)[2]
      pred.weights <- matrix(NA, M, K)
      for (m in 1:M) {
        for (k in 1:K) {
          pred.weights[m, k] <- exp(-sum((pred.coords[m, ] -
                                            mc.locations[k, ])^2)/(2 * lambda.w))
        }
        pred.weights[m, ] <- pred.weights[m, ]/sum(pred.weights[m,
                                                                ])
      }
      pred.kernel.ellipses <- array(0, dim = c(2, 2, M))
      for (m in 1:M) {
        for (k in 1:K) {
          pred.kernel.ellipses[, , m] <- pred.kernel.ellipses[,
                                                              , m] + pred.weights[m, k] * mc.kernels[, , k]
        }
      }
      if (ns.nugget == TRUE) {
        mc.nuggets <- mc.MLEs$tausq
        pred.nuggets <- rep(0, M)
        for (m in 1:M) {
          for (k in 1:K) {
            pred.nuggets[m] <- pred.nuggets[m] + pred.weights[m,
                                                              k] * mc.nuggets[k]
          }
        }
      }
      if (ns.nugget == FALSE) {
        pred.nuggets <- rep(tausq.est, M)
      }
      if (ns.variance == TRUE) {
        obs.variance <- sigmasq.est
        mc.variance <- mc.MLEs$sigmasq
        pred.variance <- rep(0, M)
        for (m in 1:M) {
          for (k in 1:K) {
            pred.variance[m] <- pred.variance[m] + pred.weights[m,
                                                                k] * mc.variance[k]
          }
        }
      }
      if (ns.variance == FALSE) {
        obs.variance <- rep(sigmasq.est, N)
        pred.variance <- rep(sigmasq.est, M)
      }
      cat("Calculating the cross-correlations. (This step can be time consuming, depending on the number of prediction locations.)")
      Scale.cross <- matrix(NA, M, N)
      Dist.cross <- matrix(NA, M, N)
      cat("\n")
      cat("Progress: ")
      for (i in 1:N) {
        Kerneli <- kernel.ellipses[, , i]
        det_i <- Kerneli[1, 1] * Kerneli[2, 2] - Kerneli[1, 2] *
          Kerneli[2, 1]
        for (j in 1:M) {
          Kernelj <- pred.kernel.ellipses[, , j]
          det_j <- Kernelj[1, 1] * Kernelj[2, 2] - Kernelj[1,
                                                           2] * Kernelj[2, 1]
          avg_ij <- 0.5 * (Kerneli + Kernelj)
          Uij <- chol(avg_ij)
          det_ij <- avg_ij[1, 1] * avg_ij[2, 2] - avg_ij[1,
                                                         2] * avg_ij[2, 1]
          vec_ij <- backsolve(Uij, (coords[i, ] - pred.coords[j,
                                                              ]), transpose = TRUE)
          Scale.cross[j, i] <- sqrt(sqrt(det_i * det_j)/det_ij)
          Dist.cross[j, i] <- sqrt(sum(vec_ij^2))
        }
        if (i%%floor(N/10) == 0) {
          cat(100 * round(i/N, 2), "% ... ", sep = "")
        }
      }
      cat("100% \n")
      Unscl.cross <- cov.spatial(Dist.cross, cov.model = cov.model,
                                 cov.pars = c(1, 1), kappa = kappa.MLE)
      NS.cross.corr <- Scale.cross * Unscl.cross
      CrossCov <- diag(sqrt(pred.variance)) %*% NS.cross.corr %*%
        diag(sqrt(obs.variance))
      pred.means <- matrix(NA, M, p)
      for (i in 1:p) {
        pred.means[, i] <- Xpred %*% beta.MLE + (CrossCov) %*%
          Cov.mat.inv %*% (object$data[, i] - Xmat %*%
                             beta.MLE)
      }
      pred.SDs <- sqrt((pred.variance + pred.nuggets) - diag((CrossCov) %*%
                                                               Cov.mat.inv %*% t(CrossCov)))
      output <- list(pred.means = pred.means, pred.SDs = pred.SDs)
      return(output)
    }
  }
}

#======================================================================================
# Fit the anisotropic model
#======================================================================================
# The Aniso_fit() function estimates the parameters of the anisotropic spatial
# model. Required inputs are the observed data and locations (a geoR object with
# $coords and $data). Optional inputs include the covariance model (exponential is
# the default) and the mean model.
#======================================================================================
#ROxygen comments ----
#' Fit the stationary spatial model
#'
#' \code{Aniso_fit} estimates the parameters of the stationary spatial model.
#' Required inputs are the observed data and locations (a geoR object
#' with $coords and $data). Optional inputs include the covariance model
#' (exponential is the default).
#'
#' @param geodata A list containing elements \code{coords} and \code{data} as
#' described next. Typically an object of the class "\code{geodata}", although
#' a geodata object only allows \code{data} to be a vector (no replicates).
#' If not provided, the arguments \code{coords} and \code{data} must be
#' provided instead.
#' @param sp.SPDF A "\code{SpatialPointsDataFrame}" object, which contains the
#' spatial coordinates and additional attribute variables corresponding to the
#' spatoal coordinates
#' @param coords An N x 2 matrix where each row has the two-dimensional
#' coordinates of the N data locations. By default, it takes the \code{coords}
#' component of the argument \code{geodata}, if provided.
#' @param data A vector or matrix with N rows, containing the data values.
#' Inputting a vector corresponds to a single replicate of data, while
#' inputting a matrix corresponds to replicates. In the case of replicates,
#' the model assumes the replicates are independent and identically
#' distributed.
#' @param cov.model A string specifying the model for the correlation
#' function; following \code{geoR}, defaults to \code{"exponential"}.
#' Options available in this package are: "\code{exponential}",
#' \code{"cauchy"}, \code{"matern"}, \code{"circular"}, \code{"cubic"},
#' \code{"gaussian"}, \code{"spherical"}, and \code{"wave"}. For further
#' details, see documentation for \code{\link[geoR]{cov.spatial}}.
#' @param mean.model An object of class \code{\link[stats]{formula}},
#' specifying the mean model to be used. Defaults to an intercept only.
#'
#' @param local.pars.LB,local.pars.UB Optional vectors of lower and upper
#' bounds, respectively, used by the \code{"L-BFGS-B"} method option in the
#' \code{\link[stats]{optim}} function for the local parameter estimation.
#' Each vector must be of length five,
#' containing values for lam1, lam2, tausq, sigmasq, and nu. Default for
#' \code{local.pars.LB} is \code{rep(1e-05,5)}; default for
#' \code{local.pars.UB} is \code{c(max.distance/2, max.distance/2, 4*resid.var, 4*resid.var, 100)},
#' where \code{max.distance} is the maximum interpoint distance of the
#' observed data and \code{resid.var} is the residual variance from using
#' \code{\link[stats]{lm}} with \code{mean.model}.
#'
#' @param local.ini.pars Optional vector of initial values used by the
#' \code{"L-BFGS-B"} method option in the \code{\link[stats]{optim}}
#' function for the local parameter estimation. The vector must be of length
#' five, containing values for lam1, lam2, tausq, sigmasq, and nu. Defaults
#' to \code{c(max.distance/10, max.distance/10, 0.1*resid.var, 0.9*resid.var, 1)},
#' where \code{max.distance} is the maximum interpoint distance of the
#' observed data and \code{resid.var} is the residual variance from using
#' \code{\link[stats]{lm}} with \code{mean.model}.
#'
#' @return A list with the following components:
#' \item{MLEs.save}{Table of local maximum likelihood estimates for each
#' mixture component location.}
#' \item{data}{Observed data values.}
#' \item{beta.GLS}{Vector of generalized least squares estimates of beta,
#' the mean coefficients.}
#' \item{beta.cov}{Covariance matrix of the generalized least squares
#' estimate of beta.}
#' \item{Mean.coefs}{"Regression table" for the mean coefficient estimates,
#' listing the estimate, standard error, and t-value.}
#' \item{Cov.mat}{Estimated covariance matrix (\code{N.obs} x \code{N.obs})
#' using all relevant parameter estimates.}
#' \item{Cov.mat.inv}{Inverse of \code{Cov.mat}, the estimated covariance
#' matrix (\code{N.obs} x \code{N.obs}).}
#' \item{aniso.pars}{Vector of MLEs for the anisotropy parameters lam1,
#' lam2, eta.}
#' \item{aniso.mat}{2 x 2 anisotropy matrix, calculated from
#' \code{aniso.pars}.}
#' \item{tausq.est}{Scalar maximum likelihood estimate of tausq (nugget
#' variance).}
#' \item{sigmasq.est}{Scalar maximum likelihood estimate of sigmasq
#' (process variance).}
#' \item{kappa.MLE}{Scalar maximum likelihood estimate for kappa (when
#' applicable).}
#' \item{cov.model}{String; the correlation model used for estimation.}
#' \item{coords}{N x 2 matrix of observation locations.}
#' \item{global.loglik}{Scalar value of the maximized likelihood from the
#' global optimization (if available).}
#' \item{Xmat}{Design matrix, obtained from using \code{\link[stats]{lm}}
#' with \code{mean.model}.}
#'
#' @examples
#' # Using iid standard Gaussian data
#' aniso.fit <- Aniso_fit( coords = cbind(runif(100), runif(100)),
#' data = rnorm(100) )
#'
#'
#' @export
#' @importFrom geoR cov.spatial
#' @importFrom StatMatch mahalanobis.dist
#' @importFrom stats lm
#' @importFrom stats optim

Aniso_fit <- function( geodata = NULL, sp.SPDF = NULL,
                       coords = geodata$coords, data = geodata$data,
                       cov.model = "exponential", mean.model = data ~ 1,
                       local.pars.LB = NULL, local.pars.UB = NULL,
                       local.ini.pars = NULL ){

  #===========================================================================
  # Formatting for coordinates/data
  #===========================================================================
  if( is.null(geodata) == FALSE ){
    if( class(geodata) != "geodata" ){
      cat("\nPlease use a geodata object for the 'geodata = ' input.\n")
    }
    coords <- geodata$coords
    data <- geodata$data
  }
  if( is.null(sp.SPDF) == FALSE ){
    if( class(sp.SPDF) != "SpatialPointsDataFrame" ){
      cat("\nPlease use a SpatialPointsDataFrame object for the 'sp.SPDF = ' input.\n")
    }
    geodata <- geoR::as.geodata( sp.SPDF )
    coords <- geodata$coords
    data <- geodata$data
  }

  N <- dim(coords)[1]
  data <- as.matrix(data)
  p <- dim(data)[2]

  # Make sure cov.model is one of the permissible options
  if( cov.model != "cauchy" & cov.model != "matern" & cov.model != "circular" &
      cov.model != "cubic" & cov.model != "gaussian" & cov.model != "exponential" &
      cov.model != "spherical" & cov.model != "wave" ){
    cat("Please specify a valid covariance model (cauchy, matern, circular, cubic, gaussian, exponential, spherical, or wave).")
  }

  #===========================================================================
  # Calculate the design matrix
  #===========================================================================
  OLS.model <- lm( mean.model, x=TRUE )

  Xmat <- matrix( unname( OLS.model$x ), nrow=N )

  #===========================================================================
  # Specify lower, upper, and initial parameter values for optim()
  #===========================================================================
  lon_min <- min(coords[,1])
  lon_max <- max(coords[,1])
  lat_min <- min(coords[,2])
  lat_max <- max(coords[,2])

  max.distance <- sqrt( sum((c(lon_min,lat_min) - c(lon_max,lat_max))^2))

  if( p > 1 ){
    ols.sigma <- NULL
    for(i in 1:length(names(summary(OLS.model)))){
      ols.sigma <- c( ols.sigma, summary(OLS.model)[[i]]$sigma )
    }
    resid.var <- (max(ols.sigma))^2
  }
  if( p == 1 ){
    resid.var <- summary(OLS.model)$sigma^2
  }

  #=================================
  # Lower limits for optim()
  #=================================
  if( is.null(local.pars.LB) == TRUE ){
    lam1.LB <- 1e-05
    lam2.LB <- 1e-05
    tausq.local.LB <- 1e-05
    sigmasq.local.LB <- 1e-05
    kappa.local.LB <- 1e-05
  }
  if( is.null(local.pars.LB) == FALSE ){
    lam1.LB <- local.pars.LB[1]
    lam2.LB <- local.pars.LB[2]
    tausq.local.LB <- local.pars.LB[3]
    sigmasq.local.LB <- local.pars.LB[4]
    kappa.local.LB <- local.pars.LB[5]
  }

  #=================================
  # Upper limits for optim()
  #=================================
  if( is.null(local.pars.UB) == TRUE ){
    lam1.UB <- max.distance/4
    lam2.UB <- max.distance/4
    tausq.local.UB <- 4*resid.var
    sigmasq.local.UB <- 4*resid.var
    kappa.local.UB <- 30
  }
  if( is.null(local.pars.UB) == FALSE ){
    lam1.UB <- local.pars.UB[1]
    lam2.UB <- local.pars.UB[2]
    tausq.local.UB <- local.pars.UB[3]
    sigmasq.local.UB <- local.pars.UB[4]
    kappa.local.UB <- local.pars.UB[5]
  }

  #=================================
  # Initial values for optim()
  #=================================

  if( is.null(local.ini.pars) == TRUE ){
    lam1.init <- max.distance/10
    lam2.init <- max.distance/10
    tausq.local.init <- 0.1*resid.var
    sigmasq.local.init <- 0.9*resid.var
    kappa.local.init <- 1
  }
  if( is.null(local.ini.pars) == FALSE ){
    lam1.init <- local.ini.pars[1]
    lam2.init <- local.ini.pars[2]
    tausq.local.init <- local.ini.pars[3]
    sigmasq.local.init <- local.ini.pars[4]
    kappa.local.init <- local.ini.pars[5]
  }

  #===========================================================================
  # MLEs
  #===========================================================================

  cat("Estimating the variance/covariance parameters. \n")

  # Covariance models with the kappa parameter
  if( cov.model == "matern" || cov.model == "cauchy" ){

    # Calculate a locally anisotropic covariance
    # Parameter order is lam1, lam2, eta, tausq, sigmasq, kappa
    anis.model.kappa <- make_aniso_loglik_kappa( locations = coords,
                                                 cov.model = cov.model,
                                                 data = data,
                                                 Xmat = Xmat )

    MLEs <- optim( c(lam1.init, lam2.init, pi/4, tausq.local.init,
                     sigmasq.local.init, kappa.local.init ),
                   anis.model.kappa,
                   method = "L-BFGS-B",
                   lower=c( lam1.LB, lam2.LB, 0, tausq.local.LB, sigmasq.local.LB ),
                   upper=c( lam1.UB, lam2.UB, pi/2,
                            tausq.local.UB, sigmasq.local.UB, kappa.local.UB ) )
    if( MLEs$convergence != 0 ){
      cat("There was an error with optim(): \n", MLEs$message, "\n")
    }
    MLE.pars <- MLEs$par
  }
  # Covariance models without the kappa parameter
  if( cov.model != "matern" & cov.model != "cauchy" ){

    # Calculate a locally anisotropic covariance
    # Parameter order is lam1, lam2, eta, tausq, sigmasq
    anis.model <- make_aniso_loglik( locations = coords,
                                     cov.model = cov.model,
                                     data = data,
                                     Xmat = Xmat )

    MLEs <- optim( c( lam1.init, lam2.init, pi/4, tausq.local.init,
                      sigmasq.local.init ),
                   anis.model,
                   method = "L-BFGS-B",
                   lower=c( lam1.LB, lam2.LB, 0, tausq.local.LB, sigmasq.local.LB ),
                   upper=c( lam1.UB, lam2.UB, pi/2,
                            tausq.local.UB, sigmasq.local.UB ) )
    if( MLEs$convergence != 0 ){
      cat("There was an error with optim(): \n", MLEs$message, "\n")
    }

    MLE.pars <- c(MLEs$par, NA)
  }

  # Save all MLEs
  # Parameter order: lam1, lam2, eta, tausq, sigmasq, mu, kappa
  MLEs.save <- MLE.pars

  # Put the MLEs into a data frame
  names(MLEs.save) <- c( "lam1", "lam2", "eta", "tausq", "sigmasq", "kappa" )
  MLEs.save <- data.frame( t(MLEs.save) )

  aniso.mat <- kernel_cov( MLE.pars[1:3] )

  global.lik <-  -MLEs$value

  #===========================================================================
  # Calculate the covariance matrix using the MLEs
  #===========================================================================
  distances <- mahalanobis.dist( data.x = coords, vc = aniso.mat )
  NS.cov <- MLE.pars[5] * cov.spatial( distances, cov.model = cov.model,
                                       cov.pars = c(1,1), kappa = MLE.pars[6] )
  Data.Cov <- NS.cov + diag(rep(MLE.pars[4], N))

  #===========================================================================
  # Calculate the GLS estimate of beta
  #===========================================================================
  Data.Cov.inv <- solve(Data.Cov)
  beta.cov <- chol2inv( chol(t(Xmat)%*%Data.Cov.inv%*%Xmat) )/p
  beta.GLS <- (p*beta.cov %*% t(Xmat) %*% Data.Cov.inv %*% data %*% rep(1,p))/p
  Mean.coefs <- data.frame( Estimate = beta.GLS,
                            Std.Error = sqrt(diag(beta.cov)),
                            t.val = beta.GLS/sqrt(diag(beta.cov)) )


  #===========================================================================
  # Output
  #===========================================================================
  output <- list( MLEs.save = MLEs.save,
                  data = data,
                  beta.GLS = beta.GLS,
                  beta.cov = beta.cov,
                  Mean.coefs = Mean.coefs,
                  Cov.mat = Data.Cov,
                  Cov.mat.inv = Data.Cov.inv,
                  aniso.pars =  MLE.pars[1:3],
                  aniso.mat = aniso.mat,
                  tausq.est = MLE.pars[4],
                  sigmasq.est = MLE.pars[5],
                  kappa.MLE = MLE.pars[6],
                  cov.model = cov.model,
                  coords = coords,
                  Xmat = Xmat,
                  global.loglik = global.lik )

  class(output) <- "Aniso"

  return(output)

}
#======================================================================================
# Calculate predictions using the output of Aniso.fit()
#======================================================================================
# Using the output from Aniso.fit(), calculate the kriging predictors
# and kriging standard errors for prediction locations of interest.
#======================================================================================
#ROxygen comments ----
#' Obtain predictions at unobserved locations for the stationary
#' spatial model.
#'
#' \code{predict.Aniso} calculates the kriging predictor and corresponding
#' standard errors at unmonitored sites.
#'
#' @param object An "Aniso" object, from \code{Aniso_fit}.
#' @param pred.coords Matrix of locations where predictions are required.
#' @param pred.covariates Matrix of covariates for the prediction locations,
#' NOT including an intercept. The number of columns for this matrix must
#' match the design matrix from \code{mean.model} in \code{\link{NSconvo_fit}}.
#' Defaults to an intercept only.
#'
#' @return A list with the following components:
#' \item{pred.means}{Vector of the kriging predictor, for each location in
#' \code{pred.coords}.}
#' \item{pred.SDs}{Vector of the kriging standard errors, for each location
#' in \code{pred.coords}.}
#' @param ... additional arguments affecting the predictions produced.
#'
#' @examples
#' \dontrun{
#' pred.S <- predict( Aniso.obj,
#' pred.coords = cbind(runif(300),runif(300)) )
#' }
#'
#' @export
#' @importFrom geoR cov.spatial
#' @importFrom StatMatch mahalanobis.dist

predict.Aniso <- function(object, pred.coords, pred.covariates = NULL,
                          ... )
{
  if( !inherits(object, "Aniso") ){
    warning("Object is not of type Aniso")
  }
  else{
    {
      M <- dim(pred.coords)[1]
      beta.MLE <- as.matrix(object$beta.GLS)
      tausq.est <- object$tausq.est
      sigmasq.est <- object$sigmasq.est
      kappa.MLE <- object$kappa.MLE
      Cov.mat.inv <- object$Cov.mat.inv
      data <- object$data
      N <- length(object$data)
      coords <- object$coords
      cov.model <- object$cov.model
      Xmat <- object$Xmat
      aniso.mat <- object$aniso.mat
      if (is.null(pred.covariates) == TRUE) {
        Xpred <- rep(1, M)
      }
      if (is.null(pred.covariates) == FALSE) {
        Xpred <- cbind(rep(1, M), pred.covariates)
      }
      data <- matrix(data, nrow = N)
      p <- dim(data)[2]
      CC.distances <- mahalanobis.dist(data.x = pred.coords, data.y = coords,
                                       vc = aniso.mat)
      CrossCov <- sigmasq.est * cov.spatial(CC.distances, cov.model = cov.model,
                                            cov.pars = c(1, 1), kappa = kappa.MLE)
      pred.means <- matrix(NA, M, p)
      for (i in 1:p) {
        pred.means[, i] <- Xpred %*% beta.MLE + ((CrossCov) %*%
                                                   Cov.mat.inv %*% (object$data[, i] - Xmat %*%
                                                                      beta.MLE))
      }
      pred.SDs <- sqrt(rep(tausq.est + sigmasq.est, M) - diag((CrossCov) %*%
                                                                Cov.mat.inv %*% t(CrossCov)))
      output <- list(pred.means = pred.means, pred.SDs = pred.SDs)
      return(output)    }
  }
}
