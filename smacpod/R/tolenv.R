# Tolerance envelopes for log ratio of spatial densities.
# 
# \code{tolenv} determines tolerance envelopes for the log ratio of spatial densities of cases and controls.  Specifically, the function identifies where the observed log ratio of spatial densities exceeds what is expected under the random labeling hypothesis.  Results can be easily plotted using the contour or image functions.
# 
# \code{alternative ="two.sided"} identifies locations where the observed log ratio of spatial densities is below and above, respectively, the (1-level)/2 and 1 - (1-level)/2 quantiles of log ratios of spatial densities simulated under the random labeling hypothesis.  "greater" finds where the observed ratio exceeds the "level" quantile.  "lower" finds where the observed ratio exceeds the 1 - level quantile.
# 
# The \code{z} argument of the list returned has a -1 for locations where the observed log ratio of spatial densities is below the tolerance envelope, a 0 for locations within the tolerance envelope, and a 1 for locations where the log ratio of spatial densities exceeds the tolerance envelope.
# 
# @param object An \code{im} object from the \code{logrr} function.
# @param level Confidence level.  Should be a number between 0 and 1.  Default is 0.95.
# @param alternative Default is "two.sided".  Can also be "greater" or "lower".  
# 
# @return Returns an \code{tolenv} object.  This is just a list with components \code{x}, \code{y}, and \code{z} the can be used with the \code{image} or \code{contour} functions quite easily.
# @author Joshua French
# @import spatstat
# @references Waller, L.A. and Gotway, C.A. (2005).  Applied Spatial Statistics for Public Health Data.  Hoboken, NJ: Wiley.
tolenv = function(object, level = 0.90, alternative = "two.sided")
{
  alpha = 1 - level
  if(alternative == "two.sided")
  {
    tol = apply(object$simr, c(1, 2), quantile, 
                prob = c(alpha/2, 1-alpha/2), na.rm = TRUE)
  }
  else if(alternative == "lower")
  {
    tol = apply(object$simr, c(1, 2), quantile, 
                prob = c(1 - level, 1), na.rm = TRUE)
  }
  else
  {
    tol = apply(object$simr, c(1, 2), quantile, 
                prob = c(0, level), na.rm = TRUE)
  }
  above = (object$simr[,,1] > tol[2,,]) + 0
  below = -1*(object$simr[,,1] < tol[1,,])
  both = above + below
  return(spatstat::im(mat = both, xcol = object$xcol, yrow = object$yrow))
}