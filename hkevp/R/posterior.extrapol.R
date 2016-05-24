#' @title
#' Spatial extrapolation with the HKEVP.
#' 
#' @author
#' Quentin Sebille
#' 
#' 
#' 
#' @description 
#' Extrapolation of the marginal GEV distribution at a set of ungauged sites (targets), given the results from the MCMC procedure \code{hkevp.fit}. See details.
#' 
#' 
#' 
#' 
#'
#' @param fit
#' A named list.
#' Output from the \code{hkevp.fit} procedure.
#' 
#' @param targets
#' A matrix of real values.
#' Coordinates of the sites where the marginal GEV distribution should be extrapolated. Each row corresponds to a target position and each column to a coordinate.
#' 
#' @param targets.covariates
#' A matrix of real values.
#' Spatial covariates associated with each target. Each row corresponds to a target and each column to a spatial covariate. By default, only the intercept and the target coordinates are used. An error occurs if number of targets covariates and sites covariates do not match.
#' 
#' @param sites.covariates
#' A matrix of real values.
#' Spatial covariates associated to the sites. Each row corresponds to a site and each column to a spatial covariate. By default, only the intercept and the target coordinates are used. An error occurs if number of targets covariates and sites covariates do not match.
#'
#'
#'
#'
#' @return
#' A named list with three elements: \code{mu}, \code{sigma}, \code{xi}, each one corresponding to a GEV parameter. Each element is a matrix where each column corresponds to a position and each row to a state of the Markov chain.
#' 
#' 
#' 
#' @details 
#' Spatial extrapolation of the GEV distribution at target positions \eqn{(s^*_1, ..., s^*_k)} is performed with a simple kriging procedure at each MCMC step on the spatial processes associated to the GEV parameters.
#' 
#' Note that if you have furnished the spatial covariates as arguments in the function \code{hkevp.fit}, they should have been standardized.
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
#' # HKEVP fit:
#' fit <- hkevp.fit(ysim, sites, knots, niter = 100, nburn = 50, quiet = FALSE)
#' 
#' ## Extrapolation:
#' targets <- matrix(1.5, 1, 2)
#' gev.targets <- posterior.extrapol(fit, targets)
#' 
#' ## True vs predicted:
#' predicted <- sapply(gev.targets, median)
#' sd.predict <- sapply(gev.targets, sd)
#' true <- c(targets[,1]*10, sigma, xi)
#' 
#' 
#' 
posterior.extrapol <- function(fit, targets, targets.covariates, sites.covariates) {
  # Default value of targets covariates and test for compatibility with sites
  if (missing(targets.covariates)) targets.covariates <- cbind(1, targets)
  if (missing(sites.covariates)) sites.covariates <- fit$spatial.covar
  if (ncol(sites.covariates) != ncol(fit$spatial.covar))
    stop("Argument sites.covariates does not match with the output of hkevp.fit!")
  if (ncol(targets.covariates) != ncol(sites.covariates))
    stop("Spatial covariates does not match between sites and targets!")
  
  all.covariates <- rbind(sites.covariates, targets.covariates)
  
  
  # Useful parameters
  nsites <- nrow(fit$sites)
  ntargets <- nrow(targets)
  all.coord <- rbind(fit$sites, targets)
  H <- as.matrix(dist(all.coord))
  GEV.vary <- fit$spatial$vary
  nstep <- fit$nstep
  
  
  # Initializing the results
  TABLE <- matrix(NA, nstep, ntargets)
  result <- list(mu = TABLE, gamma = TABLE, xi = TABLE)
  
  # Covariance function depending on the family
  if (fit$latent.correlation.type == "gauss") covariance.fun <- function(h, sill, range) sill * exp(-1/2*(h/range) ^ 2)
  if (fit$latent.correlation.type == "expo") covariance.fun <- function(h, sill, range) sill * exp(-h / range)
  if (fit$latent.correlation.type == "mat32") covariance.fun <- function(h, sill, range) sill * (1 + sqrt(3)*h/range) * exp(-sqrt(3)*h/range)
  if (fit$latent.correlation.type == "mat52") covariance.fun <- function(h, sill, range) sill * (1 + sqrt(5)*h/range  + (5/3)*(h/range) ^ 2) * exp(-sqrt(5)*h/range)
  
  
  ## Simple kriging function
  simple.kriging <- function(obs, mean, covar.mat) {
    n.obs <- length(obs)
    as.vector(mean[-(1:n.obs)] + covar.mat[-(1:n.obs),1:n.obs] %*%
                solve(covar.mat[1:n.obs,1:n.obs]) %*%
                (obs - mean[1:n.obs]))
  }
  
  
  ## Kriging the GEV parameters at each MCMC state after burn-in
  ## Scale is transformed to log for kriging, then retransformed for the returned value of the function.
  for (iter in 1:nstep) {
    
    # Covariance functions for the GEV parameters
    sills <- fit$spatial$sills[iter,]
    ranges <- fit$spatial$ranges[iter,]
    mu.cov <- covariance.fun(h = H, sill = sills[1], range = ranges[1])
    gamma.cov <- covariance.fun(h = H, sill = sills[2], range = ranges[2])
    xi.cov <- covariance.fun(h = H, sill = sills[3], range = ranges[3])
    
    # Mean of the GEV parameters
    mu.mean <- all.covariates %*% fit$spatial$beta[iter,,1]
    gamma.mean <- all.covariates %*% fit$spatial$beta[iter,,2]
    xi.mean <- all.covariates %*% fit$spatial$beta[iter,,3]
    
    # Kriging procedure
    mu.krig <- simple.kriging(obs = fit$GEV[,1,iter], mean = mu.mean, covar.mat = mu.cov)
    gamma.krig <- simple.kriging(obs = log(fit$GEV[,2,iter]), mean = gamma.mean, covar.mat = gamma.cov)
    xi.krig <- simple.kriging(obs = fit$GEV[,3,iter], mean = xi.mean, covar.mat = xi.cov)
    if (!GEV.vary[1]) mu.krig <- rep(fit$GEV[1,1,iter], ntargets)
    if (!GEV.vary[2]) gamma.krig <- rep(log(fit$GEV[1,2,iter]), ntargets)
    if (!GEV.vary[3]) xi.krig <- rep(fit$GEV[1,3,iter], ntargets)
    
    # Saving into result table
    result$mu[iter,] <- mu.krig
    result$gamma[iter,] <- gamma.krig
    result$xi[iter,] <- xi.krig
    
  }
  
  result$gamma <- exp(result$gamma)
  names(result)[2] <- "sigma"
  
  
  
  return(result)
}