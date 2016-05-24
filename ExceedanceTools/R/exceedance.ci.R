#' Return confidence region
#' 
#' \code{exceedance.ci} returns a confidence set for an exceedance region or contour line.
#' 
#' @param statistic.sim.obj An object returned from the \code{statistic.sim} function.
#' @param conf.level The desired confidence level of the confidence region.
#' @param type Whether the function should return the null region or rejection region of exceedance confidence region  Options are \code{"null"} or \code{"rejection"}.  Default is \code{"null"}.
#' 
#' @return Returns a numeric vector with the set of pixels comprising the null or rejection region related to \code{statistic.sim.obj}.
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
#' spcov <- cov.sp(coords = coords, sp.type = "exponential", 
#'  sp.par = c(1, 1.5), error.var = 1/3, finescale.var = 0, pcoords = pcoords)
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
#' # Construct null and rejection sets for two scenarios
#' n90 <- exceedance.ci(statistic.sim.obj.less, conf.level = .90, type = "null")
#' r90 <- exceedance.ci(statistic.sim.obj.greater,conf.level = .90, type = "rejection")       
#' # Plot results
#' plot(pgrid, n90, col="blue", add = FALSE, xlab = "x", ylab = "y")
#' plot(pgrid, r90, col="orange", add = TRUE)
#' legend("bottomleft", 
#'  legend = c("contains true exceedance region with 90 percent confidence", 
#'    "is contained in true exceedance region with 90 percent confidence"),
#'    col = c("blue", "orange"), lwd = 10)  

#Construct null or rejection region based on observed statistic
#and simulated statistic
exceedance.ci <- function(statistic.sim.obj, conf.level = .95, type = "null")
{
  alternative <- statistic.sim.obj$alternative
  cv <- statistic.cv(statistic.sim.obj, conf.level = conf.level)
  if(alternative == "less")
  {
    if(type == "null")
    {
      set <- which(statistic.sim.obj$statistic >= cv)
    }else
    {
      set <- which(statistic.sim.obj$statistic < cv)		
    }
  }else if(alternative == "greater")
  {
    if(type == "null")
    {
      set <- which(statistic.sim.obj$statistic <= cv)
    }else
    {
      set <- which(statistic.sim.obj$statistic > cv)		
    }
  }
  else
  {
    if(type == "null")
    {
      set <- which(statistic.sim.obj$statistic <= cv)
    }else
    {
      set <- which(statistic.sim.obj$statistic > cv)		
    }
  }
  return(set)
}
