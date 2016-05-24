#======================================================================================
#Local likelihood estimation for covariance functions with spatially-varying
#parameters: the convoSPAT() package for R
#Mark D. Risser / The Ohio State University / 2014-2015
#======================================================================================

#======================================================================================
# Simulate Functions
#======================================================================================

#======================================================================================
#Function to calculate mixture component kernels
#======================================================================================
#This function calculates spatially-varying mixture component kernels using
#generalized linear models for each of the eigenvalues (lam1 and lam2) and the
#angle of rotation (eta).
#======================================================================================
#ROxygen comments ----
#' Calculate mixture component kernel matrices.
#'
#' \code{f_mc_kernels} calculates spatially-varying mixture component kernels using
#' generalized linear models for each of the eigenvalues (lam1 and lam2) and
#' the angle of rotation (eta).
#'
#' @param y.min Lower bound for the y-coordinate axis.
#' @param y.max Upper bound for the y-coordinate axis.
#' @param x.min Lower bound for the y-coordinate axis.
#' @param x.max Upper bound for the y-coordinate axis.
#' @param N.mc Number of mixture component locations.
#' @param lam1.coef Log-linear regression coefficients for lam1; the
#' coefficients correspond to the intercept, longitude, and latitude.
#' @param lam2.coef Log-linear regression coefficients for lam2; the
#' coefficients correspond to the intercept, longitude, and latitude.
#' @param logit.eta.coef Scaled logit regression coefficients for eta; the
#' coefficients correspond to the intercept, longitude, and latitude.
#'
#' @return A list with the following components:
#' \item{mc.locations}{A \code{N.mc} x 2 matrix of the mixture component
#' locations.}
#' \item{mc.kernels}{A \code{N.mc} x 2 x 2 array of kernel matrices
#' corresponding to each of the mixture component locations.}
#'
#' @examples
#' f_mc_kernels( y.min = 0, y.max = 5, x.min = 0,
#' x.max = 5, N.mc = 3^2, lam1.coef = c(-1.3, 0.5, -0.6),
#' lam2.coef = c(-1.4, -0.1, 0.2), logit.eta.coef = c(0, -0.15, 0.15) )
#'
#'
#' @export
f_mc_kernels <- function( y.min = 0, y.max = 5, x.min = 0, x.max = 5,
                             N.mc = 3^2, lam1.coef = c(-1.3, 0.5, -0.6),
                             lam2.coef = c(-1.4, -0.1, 0.2),
                             logit.eta.coef = c(0, -0.15, 0.15) ){

  #=======================================
  # mixture component knot locations
  mc.x <- seq(from = x.min + 0.5*(x.max - x.min)/floor(sqrt(N.mc)),
                 to = x.max - 0.5*(x.max - x.min)/floor(sqrt(N.mc)),
                 length = floor(sqrt(N.mc)) )
  mc.y <- seq(from = y.min + 0.5*(y.max - y.min)/floor(sqrt(N.mc)),
                 to = y.max - 0.5*(y.max - y.min)/floor(sqrt(N.mc)),
                 length = floor(sqrt(N.mc)) )
  mc.locations <- expand.grid( mc.x, mc.y )
  mc.locations <- matrix(c(mc.locations[,1], mc.locations[,2]), ncol=2, byrow=F)

  #=======================================
  # mixture component kernels
  lam1.vec <- lam2.vec <- logit.eta.vec <- rep(0,dim(mc.locations)[1])

  lam1.vec <- exp(lam1.coef[1] + lam1.coef[2]*mc.locations[,1]
                  + lam1.coef[3]*mc.locations[,2])
  lam2.vec <- exp(lam2.coef[1] + lam2.coef[2]*mc.locations[,1]
                  + lam2.coef[3]*mc.locations[,2])
  logit.eta.vec <- logit.eta.coef[1] + logit.eta.coef[2]*mc.locations[,1]
  + logit.eta.coef[3]*mc.locations[,2]
  eta.vec <- (pi/2)*(exp(logit.eta.vec)/(1 + exp(logit.eta.vec)))

  mc.kernels <- array(NA, dim=c(2,2,dim(mc.locations)[1]))
  for(i in 1:dim(mc.locations)[1]){
    mc.kernels[,,i] <- kernel_cov( c(lam1.vec[i], lam2.vec[i], eta.vec[i]))
  }

  #=======================================
  output <- list( mc.locations = mc.locations,
                  mc.kernels = mc.kernels )
  return(output)
}

#======================================================================================
#Function to simulate data from the nonstationary model, given mixture component
#kernels
#======================================================================================
#This function requires either a mixture component kernel object, from the
#function f.mc.kernels(), or a direct specification of the mixture component
#locations and mixture component kernels. The mean can be specified to be
#nonconstant, and any of the valid covariance models for the fitting function
#can be used as well. Simulated data can be on a grid (grid = TRUE) or not
#(grid = FALSE)
#======================================================================================
#ROxygen comments ----
#' Simulate data from the nonstationary model.
#'
#' \code{NSconvo_sim} simulates data from the nonstationary model, given
#' mixture component kernel matrices. The function requires either a mixture
#' component kernel object, from the function f.mc.kernels(), or a direct
#' specification of the mixture component locations and mixture component
#' kernels.
#'
#' @param grid Logical; indicates of the simulated data should fall on a
#' grid (\code{TRUE}) or not (\code{FALSE}).
#' @param y.min Lower bound for the y-coordinate axis.
#' @param y.max Upper bound for the y-coordinate axis.
#' @param x.min Lower bound for the y-coordinate axis.
#' @param x.max Upper bound for the y-coordinate axis.
#' @param N.obs Number of simulated data values.
#' @param sim.locations Optional \code{N.obs} x 2 matrix; allows the user
#' to specify the locations of the simulated data.
#' @param mc.kernels.obj Object from the \code{\link{f_mc_kernels}} function.
#' @param mc.kernels Optional specification of mixture component kernel
#' matrices.
#' @param mc.locations Optional specification of mixture component locations.
#' @param lambda.w Scalar; tuning parameter for the weight function.
#' @param tausq Scalar; true nugget variance.
#' @param sigmasq Scalar; true process variance.
#' @param beta.coefs Vector of true regression coefficients. Length must
#' match the number of columns in \code{covariates}.
#' @param kappa Scalar; true smoothness.
#' @param covariates Matrix with \code{N.obs} rows, corresponding to
#' covariate information for each of the simualted values.
#' @param cov.model A string specifying the correlation function. Options
#' available in this package are: "\code{exponential}", \code{"cauchy"},
#' \code{"matern"}, \code{"circular"}, \code{"cubic"}, \code{"gaussian"},
#' \code{"spherical"}, and \code{"wave"}. See \code{\link[geoR]{cov.spatial}}
#' for further documentation.
#'
#' @return A list with the following components:
#' \item{sim.locations}{Matrix of locations for the simulated values.}
#' \item{mc.locations}{Mixture component locations used for the simulated
#' data.}
#' \item{mc.kernels}{Mixture component kernel matrices used for the simulated
#' data.}
#' \item{kernel.ellipses}{\code{N.obs} x 2 x 2 array, containing the kernel
#' matrices corresponding to each of the simulated values.}
#' \item{Cov.mat}{True covariance matrix (\code{N.obs} x \code{N.obs})
#' corresponding to the simulated data.}
#' \item{sim.data}{Simulated data values.}
#' \item{lambda.w}{Tuning parameter for the weight function.}
#'
#' @examples
#' \dontrun{
#' NSconvo_sim( grid = TRUE, y.min = 0, y.max = 5, x.min = 0,
#' x.max = 5, N.obs = 20^2, sim.locations = NULL, mc.kernels.obj = NULL,
#' mc.kernels = NULL, mc.locations = NULL, lambda.w = NULL,
#' tausq = 0.1, sigmasq = 1, beta.coefs = 4, kappa = NULL,
#' covariates = rep(1,N.obs), cov.model = "exponential" )
#' }
#'
#' @export
#' @importFrom stats runif

NSconvo_sim <- function( grid = TRUE, y.min = 0, y.max = 5, x.min = 0,
                         x.max = 5, N.obs = 20^2, sim.locations = NULL, mc.kernels.obj = NULL,
                         mc.kernels = NULL, mc.locations = NULL, lambda.w = NULL,
                         tausq = 0.1, sigmasq = 1, beta.coefs = 4, kappa = NULL,
                         covariates = rep(1,N.obs), cov.model = "exponential" ){

  #=======================================
  # Observed locations
  if( is.null(sim.locations) == TRUE ){
    if( grid == TRUE ){
      sim.x <- seq(from = x.min, to = x.max, length = floor(sqrt(N.obs)))
      sim.y <- seq(from = y.min, to = y.max, length = floor(sqrt(N.obs)))
      sim.locations <- expand.grid(sim.x, sim.y)
      sim.locations <- matrix(c(sim.locations[,1], sim.locations[,2]), ncol=2, byrow=F)
    }
    if( grid == FALSE ){
      sim.x <- (x.max - x.min)*runif(N.obs) + x.min
      sim.y <- (y.max - y.min)*runif(N.obs) + y.min
      sim.locations <- matrix(c(sim.x, sim.y), ncol=2, byrow=F)
      covariates = cbind(rep(1, N.obs), sim.locations)
    }
  }

  if( is.null(mc.kernels) == TRUE ){
    mc.kernels <- mc.kernels.obj$mc.kernels
  }
  if( is.null(mc.locations) == TRUE ){
    mc.locations <- mc.kernels.obj$mc.locations
  }

  #=======================================
  # Set the tuning parameter, if not specified
  if( is.null(lambda.w) == TRUE ){
    lambda.w <- ( 0.5*sqrt(sum((mc.locations[1,] - mc.locations[2,])^2)) )^2
  }

  #=======================================
  # Calculate the weights for each observed location
  K <- dim(mc.locations)[1]
  N <- dim(sim.locations)[1]

  Weights <- matrix(NA, N, K)
  for(n in 1:N){
    for(k in 1:K){
      Weights[n,k] <- exp(-sum((sim.locations[n,] - mc.locations[k,])^2)/(2*lambda.w))
    }
    # Normalize the weights
    Weights[n,] <- Weights[n,]/sum(Weights[n,])
  }

  #=======================================
  # Calculate the kernel ellipses
  kernel.ellipses <- array(0, dim=c(2,2,N))

  for(n in 1:N){
    for(k in 1:K){
      kernel.ellipses[,,n] <- kernel.ellipses[,,n] + Weights[n,k]*mc.kernels[,,k]
    }
  }

  Scale.mat <- matrix(rep(NA, N^2), nrow=N)
  Dist.mat <- matrix(rep(NA, N^2), nrow=N)
  # Calculate the elements of the observed correlations matrix.
  for(i in 1:N){

    # Diagonal elements
    Kerneli <- kernel.ellipses[,,i]
    det_i <- Kerneli[1,1]*Kerneli[2,2] - Kerneli[1,2]*Kerneli[2,1]

    Scale.mat[i,i] <- 1
    Dist.mat[i,i] <- 0

    Ui <- chol(Kerneli)

    if(i < N){
      for(j in (i+1):N){ # Off-diagonal elements

        Kernelj <- kernel.ellipses[,,j]
        det_j <- Kernelj[1,1]*Kernelj[2,2] - Kernelj[1,2]*Kernelj[2,1]

        avg_ij <- 0.5 * (Kerneli + Kernelj)
        Uij <- chol(avg_ij)
        det_ij <- avg_ij[1,1]*avg_ij[2,2] - avg_ij[1,2]*avg_ij[2,1]
        vec_ij <- backsolve(Uij, (sim.locations[i,]-sim.locations[j,]), transpose = TRUE)

        Scale.mat[i,j] <- sqrt( sqrt(det_i*det_j) / det_ij )
        Dist.mat[i,j] <- sqrt(sum(vec_ij^2))

        Scale.mat[j,i] <- Scale.mat[i,j]
        Dist.mat[j,i] <- Dist.mat[i,j]

      }
    }
  }
  Unscl.corr <- geoR::cov.spatial( Dist.mat, cov.model = cov.model,
                             cov.pars = c(1,1), kappa = kappa )
  NS.corr <- Scale.mat*Unscl.corr

  #=======================================
  # Simulate data
  Cov.mat <- sigmasq*NS.corr + diag(rep(tausq,N))
  mean.vec <- as.matrix(covariates) %*% as.matrix(beta.coefs)

  sim.data <- MASS::mvrnorm( 1, mean.vec, Cov.mat )


  output <- list( sim.locations = sim.locations,
                  mc.locations = mc.locations,
                  mc.kernels = mc.kernels,
                  kernel.ellipses = kernel.ellipses,
                  Cov.mat = Cov.mat,
                  sim.data = sim.data,
                  lambda.w = lambda.w )
  return(output)


}


#======================================================================================
# Function to calculate the kernels based on (lam1, lam2, eta)
#======================================================================================
#ROxygen comments ----
#' Calculate a kernel covariance matrix.
#'
#' \code{kernel_cov} calculates a 2 x 2 matrix based on the eigendecomposition
#' components (two eigenvalues and angle of rotation).
#'
#' @param params A vector of three parameters, corresponding to
#' (lam1, lam2, eta). The eigenvalues (lam1 and lam2) must be positive.
#'
#' @return A 2 x 2 kernel covariance matrix.
#'
#' @examples
#' kernel_cov(c(1, 2, pi/3))
#'
#' @export

kernel_cov <- function(params){
  lam1 <- params[1]
  lam2 <- params[2]
  eta <- params[3]

  Pmat <- matrix(c(cos(eta),-sin(eta),sin(eta),cos(eta)),nrow=2,byrow=T)
  Dmat <- diag(c(lam1,lam2))

  Sigma <- Pmat %*% Dmat %*% t(Pmat)
  return(Sigma)
}
