#======================================================================================
# Local likelihood estimation for covariance functions with spatially-varying
# parameters: the convoSPAT() package for R
# Mark D. Risser / The Ohio State University / 2014-2015
#======================================================================================

#======================================================================================
# Summary/Plotting Functions
#======================================================================================


#======================================================================================
# Function to print parameter estimates from NSconvo_fit()
#======================================================================================
#ROxygen comments ----
#' Summarize the nonstationary model fit.
#'
#' \code{summary.NSconvo} prints relevant output from the model fitting
#' procedure.
#'
#' @param object A "NSconvo" object, from \code{NSconvo_fit}.
#' @param ... additional arguments affecting the summary produced.
#'
#' @return Text containing the model fitting results.
#'
#' @examples
#' \dontrun{
#' summary.NSconvo( NSconvo.object )
#' }
#'
#' @export

summary.NSconvo <- function( object, ... )
{
  if( !inherits(object, "NSconvo") ){
    warning("Object is not of type NSconvo.")
  }
  else{
    cat("Locally estimated mean and variance parameters: \n")
    print(round(object$MLEs.save, 5))
    cat("\n")
    cat("Estimate of the mean coefficients: ", round(object$beta.GLS,
                                                     4), "\n")
    cat("\n")
    cat("Regression table for mean coefficients: \n")
    print(round(object$Mean.coefs, 5))
    cat("\n")
    if (is.null(object$ns.nugget) == FALSE) {
      if (object$ns.nugget == FALSE) {
        cat("Estimate of the nugget variance: ", round(object$tausq.est,
                                                       4), "\n")
      }
      if (object$ns.nugget == TRUE) {
        cat("Spatially-varying nugget variance. \n Average nugget variance: ",
            round(mean(object$tausq.est), 4), "\n")
      }
    }
    if (is.null(object$ns.nugget) == TRUE) {
      cat("Estimate of the nugget variance: ", round(object$tausq.est,
                                                     4), "\n")
    }
    if (is.null(object$ns.variance) == FALSE) {
      if (object$ns.variance == FALSE) {
        cat("Estimate of the process variance: ", round(object$sigmasq.est,
                                                        4), "\n")
      }
      if (object$ns.variance == TRUE) {
        cat("Spatially-varying process variance. \n Average process variance: ",
            round(mean(object$sigmasq.est), 4), "\n")
      }
    }
    if (is.null(object$ns.variance) == TRUE) {
      cat("Estimate of the process variance: ", round(object$sigmasq.est,
                                                      4), "\n")
    }
    cat("Estimate of the smoothness: ", round(object$kappa.MLE,
                                              4))
  }
}

#======================================================================================
# Function to print parameter estimates from Aniso_fit()
#======================================================================================
#ROxygen comments ----
#' Summarize the stationary model fit.
#'
#' \code{summary.Aniso} prints relevant output from the model fitting
#' procedure.
#'
#' @param object An "Aniso" object, from \code{Aniso_fit}.
#' @param ... additional arguments affecting the summary produced.
#'
#' @return Text containing the model fitting results.
#'
#' @examples
#' \dontrun{
#' summary.Aniso( Aniso.object )
#' }
#'
#' @export
#'

summary.Aniso <- function( object, ... )
{
  if( !inherits(object, "Aniso") ){
    warning("Object is not of type Aniso.")
  }
  else{
    cat("Globally estimated mean and variance parameters: \n")
    print(round(object$MLEs.save, 5))
    cat("\n")
    cat("Estimate of the mean coefficients: ", round(object$beta.GLS,
                                                     4), "\n")
    cat("\n")
    cat("Regression table for mean coefficients: \n")
    print(round(object$Mean.coefs, 5))
    cat("\n")
  }
}



#======================================================================================
# Function to calculate the sample size for each mixture component location
# for a particular mixture component grid and fit radius
#======================================================================================
#ROxygen comments ----
#' Calculate local sample sizes.
#'
#' \code{mc_N} calculates the number of observations (sample size) that
#' fall within a certain fit radius for each mixture component location.
#'
#' @param coords A matrix of observation locations.
#' @param mc.locations A matrix of the mixture component locations to
#' use in the model fitting.
#' @param fit.radius Scalar; defines the fitting radius for local likelihood
#' estimation.
#'
#' @return A vector \code{mc.N.fit}, which summarizes the number of
#' observation locations in \code{coords} that fall within the fit radius
#' for each mixture component location.
#'
#' @examples
#' \dontrun{
#' mc_N( coords = simdata$sim.locations, mc.locations = simdata$mc.locations,
#' fit.radius = 1 )
#' }
#'
#' @export

mc_N <- function( coords, mc.locations, fit.radius ){


  K <- dim(mc.locations)[1]
  mc.N.fit <- rep(NA,K)

  for( k in 1:K ){

    temp.locs <- coords[ abs(coords[,1]-mc.locations[k,1]) <= fit.radius
                         & (abs(coords[,2] - mc.locations[k,2]) <= fit.radius), ]

    # Isolate the data/locations to be used for calculating the local kernel
    distances <- rep(NA,dim(temp.locs)[1])

    for(i in 1:dim(temp.locs)[1]){
      distances[i] <- sqrt(sum((temp.locs[i,] - mc.locations[k,])^2))
    }

    temp.locations <- temp.locs[distances <= fit.radius,]
    n.fit <- dim(temp.locations)[1]

    mc.N.fit[k] <- n.fit

  }

  return(mc.N.fit)


}


#======================================================================================
# Evaluation criteria
#======================================================================================
#ROxygen comments ----
#' Evaluation criteria
#'
#' Calculate three evaluation criteria -- continuous rank probability score
#' (CRPS), prediction mean square deviation ratio (pMSDR), and mean squared prediction
#' error (MSPE) -- comparing hold-out data and predictions.
#'
#' @param holdout.data Observed/true data that has been held out for model
#' comparison.
#' @param pred.mean Predicted mean values corresponding to the hold-out
#' locations.
#' @param pred.SDs Predicted standard errors corresponding to the hold-out
#' locations.
#'
#' @return A list with the following components:
#' \item{CRPS}{The CRPS averaged over all hold-out locations.}
#' \item{MSPE}{The mean squared prediction error.}
#' \item{pMSDR}{The prediction mean square deviation ratio.}
#'
#' @examples
#' \dontrun{
#' evaluate_CV( holdout.data = simdata$sim.data[holdout.index],
#' pred.mean = pred.NS$pred.means, pred.SDs = pred.NS$pred.SDs )
#' }
#'
#' @export
#' @importFrom stats dnorm
#' @importFrom stats pnorm

evaluate_CV <- function( holdout.data, pred.mean, pred.SDs ){

  M <- length(holdout.data)

  # CRPS
  CRPS_out <- rep(NA,M)

  for(m in 1:M){
    sigma <- pred.SDs[m]
    zscore <- (holdout.data[m] - pred.mean[m])/sigma

    CRPS_out[m] <- sigma*( (1/sqrt(pi)) - 2*dnorm(zscore) - zscore*( 2*pnorm(zscore) - 1 ) )
  }

  # MSPE
  MSPE <- mean( (holdout.data - pred.mean)^2 )

  # pMSDR
  pMSDR <- mean( (holdout.data - pred.mean)^2/pred.SDs )


  # Output
  output <- list( CRPS = mean(CRPS_out),
                  MSPE = MSPE,
                  pMSDR = pMSDR)
  return(output)
}


#======================================================================================
# Plots from the nonstationary model
#======================================================================================


#ROxygen comments ----
#' Plot from the nonstationary model.
#'
#' This function plots either the estimated anisotropy ellipses for each
#' of the mixture component locations or the estimated correlation
#' between a reference point and all other prediction locations.
#'
#' @param x A "NSconvo" object, from NSconvo_fit().
#' @param plot.ellipses Logical; indicates whether the estimated
#' ellipses should be plotted (\code{TRUE}) or estiamted correlations
#'(\code{FALSE}).
#' @param fit.radius Scalar; defines the fit radius used for the local
#' likelihood estimation.
#' @param aniso.mat 2 x 2 matrix; contains the estimated anisotropy
#' ellipse from the stationary model (for comparison).
#' @param true.mc The true mixture component ellipses, if known.
#' @param ref.loc Vector of length 2; the reference location.
#' @param all.pred.locs A matrix of all prediction locations.
#' @param grid Logical; indicates if the \code{all.pred.locs}
#' are on a rectangular grid (\code{TRUE}) or not (\code{FALSE}).
#' @param true.col Color value for the true mixture component ellipses
#' (if plotted).
#' @param aniso.col Color value for the anisotropy ellipse (if
#' plotted).
#' @param ns.col Color value for the mixture component ellipses.
#' @param plot.mc.locs Logical; indicates whether the mixture
#' component locations should be plotted (\code{TRUE}) or not
#' (\code{FALSE}).
#' @param ... Other options passed to \code{plot}.
#'
#' @return A plot of either the estimated ellipses or estimated
#' correlation is printed.
#'
#' @examples
#' \dontrun{
#' plot.NSconvo( NSconvo.object )
#' }
#'
#' @export
#' @importFrom ellipse ellipse
#' @importFrom plotrix draw.circle
#' @importFrom graphics lines
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom fields image.plot
#' @importFrom fields quilt.plot

plot.NSconvo <- function(x, plot.ellipses = TRUE,
                         fit.radius = NULL, aniso.mat = NULL, true.mc = NULL,
                         ref.loc = NULL, all.pred.locs = NULL, grid = TRUE,
                         true.col = 1, aniso.col = 4,
                         ns.col = 2, plot.mc.locs = TRUE, ... )
{
if( !inherits(x, "NSconvo") ){
  warning("Object is not of type NSconvo.")
}
else{
  if( plot.ellipses == TRUE ){

    mc.locations <- x$mc.locations
    mc.kernels <- x$mc.kernels

    K <- dim(mc.locations)[1]

    plot(ellipse(mc.kernels[, , 1], centre = mc.locations[1,], level = 0.5),
         type = "l", col = ns.col, ... )

    if( plot.mc.locs ){ points(mc.locations[1,1], mc.locations[1,2], cex = 1, pch="+" ) }

    if (is.null(aniso.mat) == FALSE) {
      lines(ellipse(aniso.mat, centre = mc.locations[1, ],
                    level = 0.5), lty = "dashed", col = aniso.col )
    }
    if (is.null(true.mc) == FALSE) {
      lines(ellipse(true.mc[,,1], centre = mc.locations[1, ],
                    level = 0.5), col = true.col )
    }
    plotrix::draw.circle(mc.locations[1, 1], mc.locations[1, 2], fit.radius,
                         lty = "dashed" )
    for (k in 2:K) {

      lines(ellipse(mc.kernels[, , k], centre = mc.locations[k,], level = 0.5), col = ns.col )

      plotrix::draw.circle(mc.locations[k, 1], mc.locations[k,2], fit.radius, lty = "dashed")

      if( plot.mc.locs ){ points(mc.locations[k,1], mc.locations[k,2], cex = 1, pch="+" ) }
      if (is.null(aniso.mat) == FALSE) {
        lines(ellipse(aniso.mat, centre = mc.locations[k,], level = 0.5), lty = "dashed",
              col = aniso.col)
      }
      if (is.null(true.mc) == FALSE) {
        lines(ellipse(true.mc[,,k], centre = mc.locations[k, ],
                      level = 0.5), col = true.col)
      }
    }
  }

  if( plot.ellipses == FALSE ){

    M <- dim(all.pred.locs)[1]

    mc.kern <- x$mc.kernels
    mc.loc <- x$mc.locations
    K <- dim(mc.loc)[1]
    lambda.w <- x$lambda.w
    kappa <- x$kappa.MLE
    cov.model <- x$cov.model
    all.pred.weights <- matrix(NA, M, K)
    for (m in 1:M) {
      for (k in 1:K) {
        all.pred.weights[m, k] <- exp(-sum((all.pred.locs[m,] - mc.loc[k, ])^2)/(2 * lambda.w))
      }
      all.pred.weights[m, ] <- all.pred.weights[m, ]/sum(all.pred.weights[m,])
    }
    pred.weight <- rep(NA, K)
    for (k in 1:K) {
      pred.weight[k] <- exp(-sum((ref.loc - mc.loc[k, ])^2)/(2 * lambda.w))
    }
    pred.weight <- pred.weight/sum(pred.weight)
    pred.kernel.ellipses <- array(0, dim = c(2, 2, M))
    for (m in 1:M) {
      for (k in 1:K) {
        pred.kernel.ellipses[, , m] <- pred.kernel.ellipses[,,m] + all.pred.weights[m, k] * mc.kern[, ,
                                                                                                    k]
      }
    }
    pred.kernel.ellipse <- matrix(rep(0, 4), nrow = 2, ncol = 2)
    for (k in 1:K) {
      pred.kernel.ellipse <- pred.kernel.ellipse + pred.weight[k] *
        mc.kern[, , k]
    }
    Scale.cross <- rep(NA, M)
    Dist.cross <- rep(NA, M)
    Kerneli <- pred.kernel.ellipse
    det_i <- Kerneli[1, 1] * Kerneli[2, 2] - Kerneli[1, 2] *
      Kerneli[2, 1]
    Ui <- chol(Kerneli)
    for (j in 1:M) {
      Kernelj <- pred.kernel.ellipses[, , j]
      det_j <- Kernelj[1, 1] * Kernelj[2, 2] - Kernelj[1,
                                                       2] * Kernelj[2, 1]
      avg_ij <- 0.5 * (Kerneli + Kernelj)
      Uij <- chol(avg_ij)
      det_ij <- avg_ij[1, 1] * avg_ij[2, 2] - avg_ij[1,
                                                     2] * avg_ij[2, 1]
      vec_ij <- backsolve(Uij, (ref.loc - all.pred.locs[j,
                                                        ]), transpose = TRUE)
      Scale.cross[j] <- sqrt(sqrt(det_i * det_j)/det_ij)
      Dist.cross[j] <- sqrt(sum(vec_ij^2))
    }
    Unscl.cross <- geoR::cov.spatial(Dist.cross, cov.model = cov.model,
                                     cov.pars = c(1, 1), kappa = kappa)
    Cov <- matrix(Scale.cross * Unscl.cross, ncol = 1)

    if (grid == TRUE) {
      grid_x <- unique(all.pred.locs[, 1])
      grid_y <- unique(all.pred.locs[, 2])
      image.plot(grid_x, grid_y, matrix(Cov, length(grid_x), length(grid_y)),
                 ... )
    }
    if (grid == FALSE) {
      quilt.plot(all.pred.locs, c(Cov), ... )
    }
  }
}
}


#======================================================================================
# Plots from the stationary model
#======================================================================================


#ROxygen comments ----
#' Plot of the estimated correlations from the stationary model.
#'
#' This function plots the estimated correlation between a reference
#' point and all other prediction locations.
#'
#' @param x An "Aniso" object, from Aniso_fit().
#' @param ref.loc Vector of length 2; the reference location.
#' @param all.pred.locs A matrix of all prediction locations.
#' @param grid Logical; indicates if the \code{all.pred.locs}
#' are on a rectangular grid (\code{TRUE}) or not (\code{FALSE}).
#' @param ... Arguments passed to \code{plot} functions.
#'
#' @return A plot of either the estimated ellipses or estimated
#' correlation is printed.
#'
#' @examples
#' \dontrun{
#' plot.Aniso( Aniso.object )
#' }
#'
#' @export
#' @importFrom fields image.plot
#' @importFrom fields quilt.plot

plot.Aniso <- function(x, ref.loc = NULL, all.pred.locs = NULL,
                       grid = TRUE, ... )
{
  if( !inherits(x, "Aniso") ){
    warning("Object is not of type Aniso")
  }
  else{

    object <- x

    lam1 <- x$aniso.pars[1]
    lam2 <- x$aniso.pars[2]
    eta <- x$aniso.pars[3]
    kappa <- x$kappa.MLE
    cov.model <- x$cov.model
    Pmat <- matrix(c(cos(eta), -sin(eta), sin(eta), cos(eta)),
                   nrow = 2, byrow = T)
    Dmat <- diag(c(lam1, lam2))
    Sigma <- Pmat %*% Dmat %*% t(Pmat)
    distances <- StatMatch::mahalanobis.dist(data.x = all.pred.locs,
                                             data.y = t(ref.loc), vc = Sigma)
    Cov <- cov.spatial(distances, cov.model = cov.model,
                       cov.pars = c(1, 1), kappa = kappa)
    Cov <- matrix(Cov, ncol = 1)

    if (grid == TRUE) {
      grid_x <- unique(all.pred.locs[, 1])
      grid_y <- unique(all.pred.locs[, 2])
      image.plot(grid_x, grid_y,
                 matrix(Cov, length(grid_x), length(grid_y)), ... )
    }
    if (grid == FALSE) {
      quilt.plot(all.pred.locs, c(Cov), ... )
    }
  }
}
