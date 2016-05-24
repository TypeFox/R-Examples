#' Log ratio of spatial densities
#' 
#' \code{logrr} computes envelopes for the log ratio of spatial density functions.  The numerator in this ratio is related to the "cases" and the denominator to the "controls".
#' 
#' @param x Point pattern (object of class "ppp").
#' @param sigma  Standard deviation of isotropic Gaussian smoothing kernel for cases. Either a numerical value, or a function that computes an appropriate value of sigma.
#' @param sigmacon Standard deviation of isotropic Gaussian smoothing kernel for controls.  Default is the same as \code{sigma}.
#' @param case The position of the name of the "case" group in levels(x$marks).  The default is 2.
#' @param nsim The number of simulated data sets from which to construct the envelopes under the random labeling hypothesis.  Default is 0.
#' @param level The level for Monte Carlo test using tolerance intervals.
#' @param alternative The direction of the test.  Default is \code{"two.sided"}.  The values \code{"less"} and \code{"greater"} are also valid.
#' @param weights	Optional weights to be attached to the points. A numeric vector, numeric matrix, or an expression.
#' @param ... Additional arguments passed to spatstat::pixellate.ppp or spatstat::as.mask to determine the pixel resolution.
#' @param bwargs A list of arguments for the bandwidth function supplied to sigma, if applicable.
#' @param edge	Logical flag: if TRUE, apply edge correction.
#' @param varcov	Variance-covariance matrix of anisotropic Gaussian kernel. Incompatible with sigma.
#' @param at	String specifying whether to compute the intensity values at a grid of pixel locations (at="pixels") or only at the points of x (at="points").
#' @param leaveoneout	Logical value indicating whether to compute a leave-one-out estimator. Applicable only when at="points".
#' @param adjust	Optional. Adjustment factor for the smoothing parameter.
#' @param diggle	Logical. If TRUE, use Diggle's edge correction, which is more accurate but slower to compute than the correction described under Details.
#' @param nreport How frequently to report progress on the simulation.  Default is 50.
#' 
#' @return The function produces an object of type \code{logrrenv}.  It's components are similar to those returned by the \code{density.ppp} function from the \code{spatstat} package, with the intensity values replaces by the log ratio of spatial densities of f and g.  Includes an array \code{simr} of dimension c(nx, ny, nsim + 1), where nx and ny are the number of x and y grid points used to estimate the spatial density.  \code{simr[,,1]} is the log ratio of spatial densities for the observed data, and the remaining \code{nsim} elements in the third dimension of the array are the log ratios of spatial densities from a new ppp simulated under the random labeling hypothesis.
#' @details The \code{two.sided} alternative test assesses whether the observed ratio of log densities deviates more what is expected under the random labelling hypothesis.  When the test is significant, this suggests that the cases and controls are clustered, relative to the either.  The \code{greater} alternative test assesses whehter the cases are more clustered than the controls.  The \code{less} alternative test assesses whether the controls are more clustered than the cases.  If the estimated density of the case or control group becomes two small, this function may produce warnings due to numerical underflow.  Increasing the bandwidth (sigma) may help.  
#' 
#' @author Joshua French (and a small chunk by the authors of the )
#' @import spatstat
#' @importFrom utils flush.console
#' @export
#' @references Waller, L.A. and Gotway, C.A. (2005).  Applied Spatial Statistics for Public Health Data.  Hoboken, NJ: Wiley.  Kelsall, Julia E., and Peter J. Diggle. "Kernel estimation of relative risk." Bernoulli (1995): 3-16.  Kelsall, Julia E., and Peter J. Diggle. "Non-parametric estimation of spatial variation in relative risk." Statistics in Medicine 14.21-22 (1995): 2335-2342.
#' @examples 
#' data(grave)
#' r = logrr(grave)
#' plot(r)
#' r2 = logrr(grave, sigma = spatstat::bw.scott)
#' plot(r2)
#' rsim = logrr(grave, nsim = 9)
#' plot(rsim)

logrr = function(x, sigma = NULL, sigmacon = NULL, case = 2, nsim = 0, level = 0.90, alternative = "two.sided", ..., bwargs = list(), weights=NULL, edge=TRUE, varcov=NULL, at="pixels", leaveoneout=TRUE, adjust=1, diggle=FALSE, nreport = 50)
{
  if(!is.element("ppp", class(x))) stop("x must be a ppp object")
  if(is.null(x$marks)) stop("x must be marked as cases or controls")
  if(!is.factor(x$marks)) stop("The marks(x) must be a factor")
  nlev = length(levels(x$marks))
  if(case < 1 || case > nlev) stop("case must be an integer between 1 and length(levels(x$marks))")
  if(nsim < 0 | !is.finite(nsim)) stop("nsim must be a non-negative integer")
  if(level <= 0 | level >= 1) stop("level must be between 0 and 1.")
  if(!is.element(alternative, c("two.sided", "greater", "lower"))) stop("alternative is not valid.")
  alpha = 1 - level

  # determine bandwidth if necessary
  if(is.function(sigma)) # use user-supplied function, if given
  {
    which_bwargs <- which(names(bwargs) %in% names(formals(sigma))[-1])
    if(length(which_bwargs) > 0 )
    {
      sigma = do.call(sigma, c(list(X = x, bwargs[which_bwargs])))
    }else
    {
      sigma = do.call(sigma, list(X = x))
    }
  }
  if(is.null(sigma)) # use spatstat::bw.relrisk if nothing given
  {
    which_bwargs <- names(bwargs) %in% names(formals(spatstat::bw.relrisk))[-1]
    if(length(which_bwargs) > 0 )
    {
      sigma = do.call(spatstat::bw.relrisk, c(list(X = x, bwargs[which_bwargs])))
    }else
    {
      sigma = do.call(spatstat::bw.relrisk, list(X = x))
    }
  }
  if(is.null(sigmacon)) sigmacon = sigma # determine sigmacon, if NULL
  
  cases = which(x$marks == levels(x$marks)[case])
  N1 = length(cases)
  r = spdensity(x = x[cases,], sigma = sigma, ..., weights = weights,
                  edge = edge, varcov = varcov, at = at, leaveoneout = leaveoneout,
                  adjust = adjust, diggle = diggle)
  g = spdensity(x = x[-cases,], sigma = sigmacon, ..., weights = weights,
                  edge = edge, varcov = varcov, at = at, leaveoneout = leaveoneout,
                  adjust = adjust, diggle = diggle)
  r$v = log(r$v) - log(g$v)
  r$tolenv = NULL
  
  if(nsim > 0)
  {
    simr = array(0, dim = c(r$dim, nsim + 1))
    simr[,,1] = r$v
    if(nreport <= nsim) cat("Simulations completed: ")
    
    for(i in 1:nsim)
    {
      cases = sample(x$n, N1)
      fsim = spdensity(x = x[cases,], sigma = sigma, ..., weights = weights,
                      edge = edge, varcov = varcov, at = at, leaveoneout = leaveoneout,
                      adjust = adjust, diggle = diggle)
      gsim = spdensity(x = x[-cases,], sigma = sigmacon, ..., weights = weights,
                      edge = edge, varcov = varcov, at = at, leaveoneout = leaveoneout,
                      adjust = adjust, diggle = diggle)
      # fsim = spdensity(x = x[cases,])
      # gsim = spdensity(x = x[-cases,])
      
      simr[,,i + 1] = log(fsim$v) - log(gsim$v)
      if((i %% nreport) == 0){ cat(paste(i,"")); utils::flush.console() }
    }
    r$simr = simr
    r$tolenv = tolenv(r, level = level, alternative = alternative)
    class(r) = c("logrrenv", class(r))
  }
  
  r$window = x$window
 
  return(r)
}


