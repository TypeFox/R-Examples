#' Kernel smoothed spatial density of point pattern
#' 
#' \code{spdensity} computes a kernel smoothed spatial density function from a point pattern.  This is essentially a slight modification of the \code{density.ppp} function from the \code{spatstat} package, which computes the spatial intensity of a point pattern.
#' 
#' @param x Point pattern (object of class "ppp").
#' @param sigma	Standard deviation of isotropic Gaussian smoothing kernel. Either a numerical value, or a function that computes an appropriate value of sigma.
#' @param weights	Optional weights to be attached to the points. A numeric vector, numeric matrix, or an expression.
#' @param ... Additional arguments passed to pixellate.ppp and as.mask to determine the pixel resolution, or passed to sigma if it is a function.
#' @param edge	Logical flag: if TRUE, apply edge correction.
#' @param varcov	Variance-covariance matrix of anisotropic Gaussian kernel. Incompatible with sigma.
#' @param at	String specifying whether to compute the intensity values at a grid of pixel locations (at="pixels") or only at the points of x (at="points").
#' @param leaveoneout	Logical value indicating whether to compute a leave-one-out estimator. Applicable only when at="points".
#' @param adjust	Optional. Adjustment factor for the smoothing parameter.
#' @param diggle	Logical. If TRUE, use Diggle's edge correction, which is more accurate but slower to compute than the correction described under Details.
#' 
#' @return This function produces an object of class \code{im} from the \code{spatstat} package, in nearly the exact same way as \code{spatstat::density.ppp}.  The difference is that the values are scaled so that a true spatial density function is produced (i.e., the function integrates to 1).  
#' @author Joshua French
#' @import spatstat
#' @export
#' @seealso \code{\link[spatstat]{density.ppp}}
#' @references Waller, L.A. and Gotway, C.A. (2005).  Applied Spatial Statistics for Public Health Data.  Hoboken, NJ: Wiley.
#' @examples 
#' data(grave)
#' contour(spdensity(grave))

spdensity = function(x, sigma = NULL, ..., weights=NULL, edge=TRUE, varcov=NULL, at="pixels", 
                       leaveoneout=TRUE, adjust=1, diggle=FALSE)
{
  d = spatstat::density.ppp(x = x, sigma = sigma, ..., weights = weights,
              edge = edge, varcov = varcov, at = at, leaveoneout = leaveoneout,
              adjust = adjust, diggle = diggle)
  d$const = spatstat::integral.im(d)
  d$v <- d$v/d$const
  class(d) <- c(class(d), "spdensity")
  return(d)
}


