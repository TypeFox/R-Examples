#' Mean detection probability.
#'
#' Calculates the mean detection probability over the range of observation.
#' @param theta The detection function parameters. A single value for halfnormal, or a vector of two values for hazard rate.
#' @param w The maximum range of observation. Objects at distance greater than w from the observer are assumed to never be recorded.
#' @param detection.function The detection function. Only "halfnormal" and "hazard" (hazard rate) are supported at present.
#' @return The mean detection probability over the distance range [0,w].
#' @references
#' Buckland S, Anderson D, Burnham K, Laake J and Borchers D (2001). Introduction to Distance Sampling: Estimating Abundance of Biological Populations. Oxford: Oxford University Press.
#'
#' Clark, R. G. (2016), "Statistical efficiency in distance sampling," PLoS One, forthcoming, www.plosone.org
#' @examples
#' calculate.mean.detection.prob(detection.function="hazard",theta=c(0.448,2),w=1) # should be 0.6

calculate.mean.detection.prob <- function(theta=theta,w,detection.function){
  integrate(f=detection.prob,lower=0,upper=w,theta=theta,rel.tol=1e-10,detection.function=detection.function)$value / w
}
