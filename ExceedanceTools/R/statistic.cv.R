#' Return critical value of distribution.
#' 
#' \code{statistic.cv} returns the critical value of the distribution of the test statistics from \code{statistic.sim} based on the specified confidence level. However, it is not recommended for general usage.  It is recommedned that the \code{exceedance.ci} function be used to automatically create confidence regions.
#' 
#' @param statistic.sim.obj An object returned from the \code{statistic.sim} function.
#' @param conf.level The desired confidence level of the confidence interval we want to construct.
#' 
#' @return Returns the desired critical value. 
#' @author Joshua French
#' @export
#' @examples
#' library(SpatialTools)
#' 
#' # Example for exceedance regions
#' 
#' set.seed(10)
#' # Load data
#' data(sdata)
#' # Create prediction grid
#' pgrid <- create.pgrid(0, 1, 0, 1, nx = 26, ny = 26)
#' pcoords <- pgrid$pgrid
#' # Create design matrices
#' coords = cbind(sdata$x1, sdata$x2)
#' X <- cbind(1, coords)
#' Xp <- cbind(1, pcoords)
#' 
#' # Generate covariance matrices V, Vp, Vop using appropriate parameters for 
#' # observed data and responses to be predicted
#' spcov <- cov.sp(coords = coords, sp.type = "exponential", sp.par = c(1, 1.5),
#'  error.var = 1/3, finescale.var = 0, pcoords = pcoords)
#' 
#' # Predict responses at pgrid locations
#' krige.obj <- krige.uk(y = as.vector(sdata$y), V = spcov$V, Vp = spcov$Vp, 
#'  Vop = spcov$Vop, X = X, Xp = Xp, nsim = 100, 
#'  Ve.diag = rep(1/3, length(sdata$y)) , method = "chol")
#'                 
#' # Simulate distribution of test statistic for different alternatives
#' statistic.sim.obj.less <- statistic.sim(krige.obj = krige.obj, level = 5, 
#'  alternative = "less")
#' statistic.sim.obj.greater <- statistic.sim(krige.obj = krige.obj, level = 5,
#'  alternative = "greater")
#' # Calculate quantiles of distribution of statistic
#' q90.less <- statistic.cv(statistic.sim.obj.less, conf.level = .90)
#' q90.greater <- statistic.cv(statistic.sim.obj.greater, conf.level = .90)
statistic.cv <- function(statistic.sim.obj, conf.level = .95)
{
  n <- length(statistic.sim.obj$statistic.sim)
  alternative <- statistic.sim.obj$alternative
  
  #Determine whether we can take the ith element of the sorted
  #statistics as the quantile, or if we need to average elements.  
  #No average is type == 1, average is type == 2.
  if((n * conf.level) == floor(n * conf.level))
  {
    type <- 1
  }else
  {
    type <- 2
  }
  
  #Determine critical value for statistic
  if(alternative == "less")
  {
    cv <- quantile(statistic.sim.obj$statistic.sim, 
                   prob = 1 - conf.level, type = type)
  }else if(alternative == "greater")
  {
    cv <- quantile(statistic.sim.obj$statistic.sim, 
                   prob = conf.level, type = type)
  }else
  {
    cv <- quantile(statistic.sim.obj$statistic.sim, 
                   prob = conf.level, type = type)
  }
  
  return(cv)
}
