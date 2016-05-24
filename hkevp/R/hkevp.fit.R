#' @title
#' Fitting procedure of the HKEVP
#' 
#' @author
#' Quentin Sebille
#' 
#' 
#' 
#' @description 
#' Metropolis-within-Gibbs algorithm that returns posterior distribution (as Markov chains) for the marginal and dependence parameters of the HKEVP.
#' The implementation has been done such that user has to provide the minimum input parameters as possible, but convergence of the chains should be assessed, using for instance \code{\link{mcmc.plot}}. Initial states, prior hyperparameters and magnitude of random walks can be set by the experimented user.
#' 
#' 
#' 
#' 
#' 
#' 
#' @param y
#' A numerical matrix of real values.
#' Observed block maxima process: each column corresponds to a site position whose coordinates are given by \code{sites} and each row represents a block of observations (e.g. a specific year).
#' 
#' @param sites
#' A numerical matrix of real values.
#' Coordinates of the sites: each row corresponds to a position and each column to a spatial coordinate.
#' 
#' @param knots
#' A numerical matrix of real values.
#' Coordinates of the knots: each row corresponds to a knot position and each column to a spatial coordinate. If missing, \code{knots} is equal to \code{sites}. See details for how to place the knots.
#' 
#' @param niter
#' A positive integer.
#' Number of iterations before the algorithm stops.
#' 
#' @param nburn
#' A positive integer.
#' Number of iterations to be discarded to assess convergence. Half of \code{niter} is taken by default.
#' 
#' @param nthin
#' A positive integer.
#' Size of the thinning, 1 by default.
#' 
#' @param trace
#' A positive integer.
#' If \code{quiet} is FALSE, the log-likelihood of the model is displayed each block of \code{trace} iterations, which serves to observe the progression of the algorithm.
#' 
#' @param spatial.covar
#' A numerical matrix.
#' Spatial covariates associated to the sites. Each row corresponds to a position and each column to a spatial covariate. The first column should be filled with ones to represent the intercept. See details.
#' 
#' @param quiet
#' Logical.
#' Whether or not the progression of the routine may be displayed. If TRUE, nothing is displayed (quiet mode) and if FALSE (default), log-likelihood is given every block of \code{trace} iterations and acceptance rates are shown at the end of the algorithm.
#' 
#' @param GEV.vary
#' A logical vector of size three.
#' Each element corresponds to a GEV parameter (respectively the location, the scale and the shape) and indicates if it is spatially-varying (TRUE) or spatially-constant (FALSE, by default for the shape parameter). In the latter case, the value of the GEV parameter is the same for each site.
#' 
#' @param latent.correlation.type
#' A character string.
#' Indicates the form of the correlation functions of the latent Gaussian processes associated to the GEV parameters. Must be one of \code{"expo"}, \code{"gauss"}, \code{"mat32"} (By default) and \code{"mat52"}, respectively corresponding to the exponential, Gaussian, Matern-3/2 and Matern-5/2 correlation functions.
#' 
#' @param GEV.init
#' A matrix of numerical values.
#' Initial values for the GEV parameters. Must have three columns and same number of rows as in \code{sites}. If missing, default values are used (see details).
#' 
#' @param alpha.init
#' A numerical value in (0,1].
#' Initial value of the dependence parameter \eqn{\alpha}.
#' 
#' @param tau.init
#' A positive numerical value.
#' Initial value of the bandwith parameter \eqn{\tau}. By default, the minimum distance between knots is taken. Should be lower than twice the maximum distance between sites (see details).
#' 
#' @param A.init
#' A positive numerical value.
#' Initial value of the random effect \eqn{A}. The same initial state is considered for each knot and for each observation of the process.
#' 
#' @param B.init
#' A numerical value in (0,1).
#' Initial value of the auxiliary variable \eqn{B} used to update the random effect \eqn{A}. See appendices in \cite{Reich and Shaby (2012)} or the method of \cite{Stephenson (2009)}.
#' 
#' @param range.init
#' A positive numerical value.
#' Initial value of the range parameter of the latent Gaussian processes associated to the GEV parameters. Must be lower than twice the maximum distance between sites. By default, one third of this maximum distance is taken as initial value. See details.
#' 
#' @param sill.init
#' A positive numerical value.
#' Initial value of the sill parameter of the latent Gaussian processes associated to the GEV parameters.
#' 
#' @param constant.GEV.prior
#' A numerical matrix with two lines and three columns.
#' Normal prior hyperparameters for the spatially-constant GEV parameters. The first (resp. second) line corresponds to the mean (resp. standard deviation). Each column corresponds to a GEV parameter.
#' 
#' @param latent.beta.prior
#' A numerical value.
#' Standard deviation of the normal prior associated to the \eqn{\beta} parameters in the mean of the latent processes.
#' 
#' @param alpha.prior
#' A two-dimensional vector of positive values.
#' Parameters of the Beta prior for \eqn{\alpha}. By default, the uniform distribution is considered.
#' 
#' @param tau.prior
#' A two-dimensional vector of positive values.
#' Parameters of the Beta prior for \eqn{\tau} over the interval \eqn{(0, 2D_{max})}, where \eqn{D_{max}} is the maximum distance between sites. By default, a Beta(2,5) distribution is considered.
#' 
#' @param sill.prior
#' A two-dimensional vector of positive values.
#' Parameters of the Inverse-Gamma conjugate prior for the sill parameter of the latent Gaussian processes. The prior is the same for all spatially-varying GEV parameters.
#' 
#' @param range.prior
#' A two-dimensional vector of positive values.
#' Parameters of the Beta prior for the range parameter of the latent Gaussian processes. The Beta prior is considered over the interval \eqn{(0, 2D_{max})}, where \eqn{D_{max}} is the maximum distance between sites (see details). The prior is the same for all spatially-varying GEV parameters.
#' 
#' @param GEV.random.walk
#' A three-dimensional vector of positive values.
#' Standard deviations of the normal random walk for each GEV parameters.
#' 
#' @param range.random.walk
#' A three-dimensional vector of positive values.
#' Standard deviations of the normal random walk for the range parameter of the latent Gaussian processes.
#' 
#' @param tau.random.walk
#' A positive value.
#' Standard deviation of the normal random walk for the bandwith parameter \eqn{\tau}.
#' 
#' @param alpha.random.walk
#' A positive value.
#' Standard deviation of the normal random walk for the dependence parameter \eqn{\alpha}.
#' 
#' @param A.random.walk
#' A positive value.
#' Standard deviation of the normal random walk for the random effect \eqn{A}.
#' 
#' @param B.random.walk
#' A positive value.
#' Standard deviation of the normal random walk for the auxiliary variable \eqn{B}.
#'
#'
#'
#'
#'
#' @details 
#' \itemize{
#' \item{The positive integer arguments \code{niter}, \code{nburn} and \code{nthin} control the length of the resulting chains in the MCMC algorithm. The routine operates over a loop of size \code{niter} and for each step, the elements of the Markov chains are updated \code{nthin} times to account for thinning procedure. At the end, the first \code{nburn} iterations are then discarded. 
#' As a result, the length of the Markov chains given by the \code{hkevp.fit} procedure is \code{niter-nburn}. The \code{hkevp.fit} function updates the Markov chains \code{niter*nthin} times.}
#' 
#' \item{The spatial covariates for each position are given by the \code{spatial.covar} matrix argument. These values are used in the linear modelling of the mean of the latent Gaussian processes associated to the GEV parameters. If this argument is given by the user, the first column must be filled with ones to represent the intercept and other columns should be standardized. If not given, the default spatial covariates used are the intercept and the sites positions.}
#' 
#' \item{Perhaps the most restrictive feature when fitting the HKEVP to observed data comes from the choice of the knots.
#' A trade-off may be done here, since too few knots may result in a poor fit while too many knots increases drastically the computational burden of the algorithm.
#' A possible advice for the choice of \code{knots} is to take a gridded network of positions over the studied region, by using for instance the routine \code{expand.grid} from R base package.
#' If this argument is not furnished, the knots will coincide with the site positions.}
#' 
#' \item{If the initial values of the GEV parameters are not given through \code{GEV.init}, the median of the observed data are taken for the location parameter (by site if \code{GEV.vary[1]} is TRUE, or over all data if \code{GEV.vary[1]} is FALSE). For the scale and shape parameters, the default initial values are respectively 1 and 0.01.}
#' 
#' \item{A change has been made from the original MCMC procedure described in \cite{Reich and Shaby (2012)} concerning the bandwidth parameter \eqn{\tau}. Instead of being defined over the positive real values and with a vague inverse-Gamma prior, we chose to restrict this parameter to the interval \eqn{(0,2D_{max}]}, where \eqn{D_{max}} is the maximum distance between sites. A Beta(2,5) prior is set by default on this interval. For this reason, the initial value \code{tau.init} should be between 0 and \eqn{2D_{max}}. A similar remark can be done for the range parameters of the latent Gaussian processes and the associated initial values \code{range.init}}
#' 
#' \item{In the appendix of \cite{Reich and Shaby (2012)} about the MCMC algorithm, an adaptive procedure taking place in the burn-in period is presented, but it has not been implemented in this function. Instead, acceptance rates can be displayed with the function \code{\link{mcmc.accept}}.}}
#'
#'
#'
#'
#' @return
#' A named list with twelve elements:
#' \itemize{
#' \item{\code{GEV}: a three-dimensional array of real values. Markov chains associated to the GEV parameters per site. The dimensions of the array correspond respectively to the sites positions, the three GEV parameters and the states of the Markov chains.}
#' \item{\code{alpha}: a column vector. Markov chain associated to the dependence parameter \eqn{\alpha}.}
#' \item{\code{tau}: a column vector. Markov chain associated to the dependence parameter \eqn{\tau}.}
#' \item{\code{A}: a three-dimensional array of real values. Markov chains associated to the positive stable random effect per site and per block. The dimensions correspond respectively to the indices of blocks, the knots positions and the states of the Markov chains.}
#' \item{\code{llik}: a column vector. Log-likelihood of the model for each step of the algorithm.}
#' \item{\code{time}: a positive value. Time (in sec) spent for the fit.}
#' \item{\code{spatial}: a named list with four elements linked to the GEV spatially-varying parameters:
#' \itemize{
#' \item{\code{vary}: a column vector of length three. Numerical version of the argument \code{GEV.vary}.}
#' \item{\code{beta}: a three-dimensional array of real values. The dimensions correspond respectively to the states of the Markov chains, the associated spatial covariates and the three GEV parameters}
#' \item{\code{sills}: a three column matrix of real values. Markov chains associated to the sills in the correlation functions of the latent Gaussian processes. Each column corresponds to one GEV parameter.}
#' \item{\code{ranges}: a three column matrix of real values. Markov chains associated to the ranges in the correlation functions of the latent Gaussian processes. Each column corresponds to one GEV parameter.}
#' }}
#' \item{\code{knots}: the set of knots.}
#' \item{\code{spatial.covar}}: the spatial covariates.
#' \item{\code{latent.correlation.type}: the type of correlation function that describe the latent Gaussian processes.}
#' \item{\code{nstep}: a positive integer. Number of states \code{niter-nburn} obtained at the end of the routine.}
#' }
#' 
#' 
#' 
#' 
#' @export
#'
#'
#' @references 
#' Reich, B. J., & Shaby, B. A. (2012). A hierarchical max-stable spatial model for extreme precipitation. The annals of applied statistics, 6(4), 1430. <DOI:10.1214/12-AOAS591>
#' 
#' Stephenson, A. G. (2009) High-dimensional parametric modelling of multivariate extreme events. Aust. N. Z. J Stat, 51, 77-88. <DOI:10.1111/j.1467-842X.2008.00528.x>
#' 
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
#' 
#' 
#' 
#' 
hkevp.fit <- function(y, sites, knots, niter, nburn, nthin, trace, spatial.covar, quiet, GEV.vary, latent.correlation.type, GEV.init, alpha.init, tau.init, A.init, B.init, range.init, sill.init, constant.GEV.prior, latent.beta.prior, alpha.prior, tau.prior, sill.prior, range.prior, GEV.random.walk, range.random.walk, tau.random.walk, alpha.random.walk, A.random.walk, B.random.walk) {
  # Default values for missing arguments and errors
  if (!is.matrix(y)) stop("y must be a matrix!")
  if (sum(colSums(!is.na(y)) == 0) > 0) stop("y must have at least one observation per column!")
  if (!is.matrix(sites)) stop("sites must be a matrix!")
  if (ncol(y) != nrow(sites)) stop("y and sites do not match!")
  if (missing(knots)) knots <- sites
  if (!is.matrix(knots)) stop("knots must be a matrix!")
  if (missing(nburn)) nburn <- floor(niter/2)
  if (nburn < 0) nburn <- 0
  if (nburn >= niter) stop("nburn must be lower than niter!")
  if (missing(nthin)) nthin <- 1
  if (nthin <= 0) nthin <- 1
  if (missing(trace)) trace <- max(floor(niter/20), 1)
  if (trace <= 0) trace <- 1
  trace <- floor(trace) # Protection against possible seg.fault!
  if (missing(spatial.covar)) spatial.covar <- cbind(1, sites)
  if (nrow(spatial.covar) != nrow(sites)) stop("sites and spatial.covar do not match!")
  if (missing(quiet)) quiet <- FALSE
  if (missing(GEV.vary)) GEV.vary <- c(TRUE, TRUE, FALSE)
  if (missing(latent.correlation.type)) latent.correlation.type <- "mat32"  # Default correlation for latent process is Matern-3/2
  if (!(latent.correlation.type %in% c("expo", "gauss", "mat32", "mat52"))) stop("Bad argument for latent.correlation.type!")
  if (missing(alpha.init)) alpha.init <- .25
  if (missing(tau.init)) tau.init <- min(dist(knots))
  if (tau.init <=0 | tau.init > max(dist(sites))*2) stop("tau.init must be between 0 and 2*Dmax!")
  if (missing(A.init)) A.init <- exp(2)
  if (missing(B.init)) B.init <- .5
  if (missing(range.init)) range.init <- max(dist(sites))/3
  if (range.init <=0 | range.init > max(dist(sites))*2) stop("range.init must be between 0 and 2*Dmax!")
  range.init <- min(max(range.init, 0), max(dist(sites))*2)
  if (missing(sill.init)) sill.init <- 1
  if (sill.init <= 0) stop("sill.init must be positive!")
  if (missing(constant.GEV.prior)) constant.GEV.prior <- rbind(rep(0,3), c(10, 1, .25))
  if (missing(latent.beta.prior)) latent.beta.prior <- 100
  if (missing(alpha.prior)) alpha.prior <- c(1,1)
  if (missing(tau.prior)) tau.prior <- c(2,5)
  if (missing(sill.prior)) sill.prior <- c(.1,.1)
  if (missing(range.prior)) range.prior <- c(2,5)
  sill.prior <- matrix(sill.prior, 2, 3)
  range.prior <- matrix(range.prior, 2, 3)
  if (missing(GEV.random.walk)) GEV.random.walk <- c(1, .1, .01)
  if (missing(range.random.walk)) range.random.walk <- rep(.02, 3)
  if (missing(tau.random.walk)) tau.random.walk <- .02
  if (missing(alpha.random.walk)) alpha.random.walk <- .01
  if (missing(A.random.walk)) A.random.walk <- 1
  if (missing(B.random.walk)) B.random.walk <- .25
  
  # Initial state of GEV parameters
  if (missing(GEV.init)) {
    GEV.init <- matrix(NA, nrow(sites), 3)
    GEV.init[,1] <- ifelse(GEV.vary[1], apply(y, 2, median, na.rm = TRUE), median(y, na.rm = TRUE))
    GEV.init[,2] <- 1
    GEV.init[,3] <- 0.01
  } else {
    if (!is.matrix(GEV.init)) stop("GEV.init must be a matrix!")
    if (any(nrow(GEV.init) != nrow(sites) & ncol(GEV.init) != 3)) stop("Bad dimensions for GEV.init, must be nx3!")
    if (sum(GEV.init[,2] <= 0) > 0) stop("Initial scale must be positive !")
  }
  GEV.init[,2] <- log(GEV.init[,2]) ## GEV scale transformed to GEV log-scale
  
  # Distance and NA matrices
  nsites <- nrow(sites)
  distances <- as.matrix(dist(rbind(sites, knots)))
  dss <- distances[1:nsites, 1:nsites]
  dsk <- distances[1:nsites, -(1:nsites)]
  na.mat <- 1 - is.na(y)  # Matrix of missing values
  
  # Missing values are replaced by 0 but ignored in the MCMC routine
  y0 <- y
  y0[is.na(y)] <- 0
  
  # Calling C++ function
  result <-  .Call('hkevp_MCMC', PACKAGE = 'hkevp', y0, sites, spatial.covar, knots, dsk, dss, niter, nburn, nthin, trace, GEV.vary, GEV.init, alpha.init, tau.init, A.init, B.init, sill.init, range.init, constant.GEV.prior, alpha.prior, tau.prior, latent.beta.prior, sill.prior, range.prior, GEV.random.walk, range.random.walk, tau.random.walk, alpha.random.walk, A.random.walk, B.random.walk, quiet, latent.correlation.type, na.mat)
  
  result$GEV[,2,] <- exp(result$GEV[,2,]) ## GEV log scale transformed to GEV scale
  result$sites <- sites
  result$knots <- knots
  result$spatial.covar <- spatial.covar
  result$latent.correlation.type <- latent.correlation.type
  result$nstep <- niter - nburn
  
  return(result)
}
