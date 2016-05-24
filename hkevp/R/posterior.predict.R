#' @title
#' Prediction of the max-stable process at target positions.
#' 
#' @author
#' Quentin Sebille
#' 
#' @description 
#' Returns the predictive posterior distribution of the block maxima process \eqn{Y(\cdot)} at a set of target sites \eqn{(s_1^*, ..., s_k^*)}, given the observations at \eqn{(s_1, ..., s_n)} and the output from \code{hkevp.fit}.
#' 
#' Two types of conditional prediction are available, as described in \cite{Shaby and Reich (2012)}. See details.
#' 
#' 
#' 
#'
#' @inheritParams posterior.extrapol
#' 
#' @param predict.type 
#' Character string.
#' Specifies the type of prediction. Must be one of "\code{kriging}" (default) or "\code{climat}". See details.
#' 
#' 
#'
#' @details 
#' The spatial prediction of \eqn{Y_t(s^*)} for a target site \eqn{s^*} and a realisation \eqn{t} of the process is described in \cite{Shaby and Reich (2012)}. This method involves a three-step procedure:
#' \enumerate{
#' \item Computation of the residual dependence process \eqn{\theta(\cdot)} at the target positions.
#' \item Computation of the conditional GEV parameters \eqn{(\mu^*,\sigma^*,\xi^*)} at the target sites. See the definition of the HKEVP in \cite{Reich and Shaby (2012)}.
#' \item Generation of \eqn{Y_t(s^*)} from an independent GEV distribution with parameters \eqn{(\mu^*,\sigma^*,\xi^*)}.
#' }
#' 
#' As sketched in \cite{Shaby and Reich (2012)}, two types of prediction are possible: the kriging-type and the climatological-type. These two types differ when the residual dependence process \eqn{\theta} is computed (first step of the prediction):
#' \itemize{
#' \item The kriging-type takes the actual value of \eqn{A} in the MCMC algorithm to compute the residual dependence process. The prediction will be the distribution of the maximum recorded at the specified targets.
#' \item The climatological-type generates \eqn{A} by sampling from the positive stable distribution with characteristic exponent \eqn{\alpha}, where \eqn{\alpha} is the actual value of the MCMC step. The prediction in climatological-type will be the distribution of what could happen in the conditions of the HKEVP dependence structure.
#' }
#'
#' Posterior distribution for each realisation \eqn{t} of the process and each target position \eqn{s^*} is represented with a sample where each element corresponds to a step of the MCMC procedure.
#'
#'
#'
#' @return
#' A three-dimensional array where:
#' \itemize{
#' \item Each row corresponds to a different realisation of the process (a block).
#' \item Each column corresponds to a target position.
#' \item Each slice corresponds to a MCMC step.
#' }
#' 
#' 
#' @export
#' 
#' @references 
#' Reich, B. J., & Shaby, B. A. (2012). A hierarchical max-stable spatial model for extreme precipitation. The annals of applied statistics, 6(4), 1430. <DOI:10.1214/12-AOAS591>
#' 
#' Shaby, B. A., & Reich, B. J. (2012). Bayesian spatial extreme value analysis to assess the changing risk of concurrent high temperatures across large portions of European cropland. Environmetrics, 23(8), 638-648. <DOI:10.1002/env.2178>
#'
#'
#'
#' @examples
#' # Simulation of HKEVP:
#' sites <- as.matrix(expand.grid(1:3,1:3))
#' knots <- sites
#' mu <- sites[,1]*10
#' sigma <- 3
#' xi <- .2
#' alpha <- .4
#' tau <- 1
#' ysim <- hkevp.rand(10, sites, knots, mu, sigma, xi, alpha, tau)
#' 
#' # HKEVP fit (omitting first site, used as target):
#' fit <- hkevp.fit(ysim, sites, knots, niter = 100, nburn = 50, quiet = FALSE)
#' 
#' # Extrapolation:
#' targets <- matrix(sites[1,], 1, 2)
#' ypred <- posterior.predict(fit, targets, predict.type = "kriging")
#' 
#' # Plot of the density and the true value for 4 first realizations:
#' # Could need more sites/iterations to be more precise
#' # par(mfrow = c(2, 2))
#' # plot(density(ypred[1,1,]), main = "Predictive density for year 1")
#' # abline(v = ysim[1,1], col = 2, lwd = 2)
#' # plot(density(ypred[2,1,]), main = "Predictive density for year 2")
#' # abline(v = ysim[2,1], col = 2, lwd = 2)
#' # plot(density(ypred[3,1,]), main = "Predictive density for year 3")
#' # abline(v = ysim[3,1], col = 2, lwd = 2)
#' # plot(density(ypred[4,1,]), main = "Predictive density for year 4")
#' # abline(v = ysim[4,1], col = 2, lwd = 2)
#' 
#' 
#' 
#' 
posterior.predict <- function(fit, targets, targets.covariates, sites.covariates, predict.type = "kriging") {
  # Catching errors
  if (missing(targets.covariates)) targets.covariates <- cbind(1, targets)
  if (missing(sites.covariates)) sites.covariates <- fit$spatial.covar
  if (ncol(sites.covariates) != ncol(fit$spatial.covar))
    stop("Argument sites.covariates does not match with the output of hkevp.fit!")
  if (ncol(targets.covariates) != ncol(sites.covariates))
    stop("Spatial covariates does not match between sites and targets!")
  if (!(predict.type %in% c("kriging", "climat")))
    stop("Wrong argument for predict.type: must be 'kriging' or 'climat'!")
  
  # Useful variables
  ntargets <- nrow(targets)
  nyear <- dim(fit$A)[1]
  nstep <- fit$nstep
  nknots <- nrow(fit$knots)
  dtk <- as.matrix(dist(rbind(targets, fit$knots)))[1:ntargets, -(1:ntargets)]
  if (ntargets == 1) dtk <- matrix(dtk, nrow = 1)
  RESULT <- array(NA, dim = c(nyear, ntargets, nstep))
  GEV.extrapol <- posterior.extrapol(fit, targets, targets.covariates, sites.covariates) # Spatial extrapolation of GEV distribution at targets
  
  # Function that generates one GEV(mu, sigma, xi) simulation
  gev.rand <- function(mu, sigma, xi) {
    if (xi == 0)
      mu - sigma*log(rexp(1))
    else
      mu + sigma*(rexp(1) ^ (-xi) - 1)/xi
  }
  
  
  # PREDICTION ----------
  if (predict.type == "kriging") {
    for (i in 1:nstep) {
      # Computation of the THETA process at targets
      omega <- exp(-dtk^2/(2*fit$tau[i]^2))
      omega <- sweep(omega, MARGIN = 1, STATS = rowSums(omega), FUN = "/")
      theta.targets <- (fit$A[,,i] %*% t(omega^(1/fit$alpha[i])))^fit$alpha[i]
      
      # Computation of conditional GEV parameters
      mu.star <- matrix(GEV.extrapol$mu[i,], nyear, ntargets, byrow = TRUE) + matrix(GEV.extrapol$sigma[i,]/GEV.extrapol$xi[i,], nyear, ntargets, byrow = TRUE) * (theta.targets ^ matrix(GEV.extrapol$xi[i,], nyear, ntargets, byrow = TRUE) - 1)
      sigma.star <- fit$alpha[i] * matrix(GEV.extrapol$sigma[i,], nyear, ntargets, byrow = TRUE) * theta.targets ^ matrix(GEV.extrapol$xi[i,], nyear, ntargets, byrow = TRUE)
      xi.star <- fit$alpha[i] * matrix(GEV.extrapol$xi[i,], nyear, ntargets, byrow = TRUE)
      
      # Generation of Y
      for (t in 1:nyear) {
        for (s in 1:ntargets) {
          RESULT[t, s, i] <- gev.rand(mu.star[t,s], sigma.star[t,s], xi.star[t,s])
        }
      }
    }
  }
  
  if (predict.type == "climat") {
    for (i in 1:nstep) {
      # Generating A from the actual value of alpha
      unif.gen <- matrix(runif(nyear*nknots, 0, pi), nyear, nknots)
      expo.gen <- matrix(rexp(nyear*nknots, 1), nyear, nknots)
      A.gen <- (sin((1 - fit$alpha[i])*unif.gen) / expo.gen) ^ ((1 - fit$alpha[i])/fit$alpha[i] ) *
        sin(fit$alpha[i]*unif.gen) / (sin(unif.gen)) ^ (1/fit$alpha[i])
      
      # Computing the THETA process at targets:
      omega <- exp(-dtk^2/(2*fit$tau[i]^2))
      omega <- sweep(omega, MARGIN = 1, STATS = rowSums(omega), FUN = "/")
      theta.targets <- (A.gen %*% t(omega^(1/fit$alpha[i])))^fit$alpha[i]
      
      
      # Computation of conditional GEV parameters
      mu.star <- matrix(GEV.extrapol$mu[i,], nyear, ntargets, byrow = TRUE) + matrix(GEV.extrapol$sigma[i,]/GEV.extrapol$xi[i,], nyear, ntargets, byrow = TRUE) * (theta.targets ^ matrix(GEV.extrapol$xi[i,], nyear, ntargets, byrow = TRUE) - 1)
      sigma.star <- fit$alpha[i] * matrix(GEV.extrapol$sigma[i,], nyear, ntargets, byrow = TRUE) * theta.targets ^ matrix(GEV.extrapol$xi[i,], nyear, ntargets, byrow = TRUE)
      xi.star <- fit$alpha[i] * matrix(GEV.extrapol$xi[i,], nyear, ntargets, byrow = TRUE)
      
      # Generation of Y
      for (t in 1:nyear) {
        for (s in 1:ntargets) {
          RESULT[t, s, i] <- gev.rand(mu.star[t,s], sigma.star[t,s], xi.star[t,s])
        }
      }
    }
  }
  
  
  return(RESULT)
}
