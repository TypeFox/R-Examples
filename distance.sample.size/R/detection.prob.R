#' Calculate detection probabilities.
#'
#' Calculates the detection probability at one or more distances.
#' @param d The distance or distances of interest.
#' @param detection.function The detection function. Only "halfnormal" and "hazard" (hazard rate) are supported at present.
#' @param theta The detection function parameters. A single value for halfnormal, or a vector of two values for hazard rate.
#' @param deriv Optional numeric value specifying whether a derivative is required. If missing, the function returns the detection probabilities at distances d.
#' If deriv is equal to 1 or 2, the derivatives of the detection function with respect to theta[deriv] at d are returned.
#' Note that the halfnormal detection function has only one parameter, so setting deriv=2 and detection.function="halfnormal" will result in an error.
#' @param stop Set to T to open a browser window (for debugging purposes)
#' @return A vector of detection probabilities corresponding to the distances in d.
#' @references
#' Buckland S, Anderson D, Burnham K, Laake J and Borchers D (2001). Introduction to Distance Sampling: Estimating Abundance of Biological Populations. Oxford: Oxford University Press.
#'
#' Clark, R. G. (2016), "Statistical efficiency in distance sampling," PLoS One, forthcoming, www.plosone.org
#' @examples
#' dvalues <- seq(from=0,to=1,by=0.001)
#' dprobs <- detection.prob(d=dvalues,detection.function="hazard",theta=c(0.448,2))
#' plot(dvalues,dprobs,type="l",ylim=c(0,1))

detection.prob <- function(d,detection.function=c("halfnormal","hazard"),theta,deriv,stop=F){
  if(stop) browser()
  # if deriv is supplied, then the derivative with respect to theta[1] is returned (if deriv=1) or the derivative with respect to
  #   theta[2] is returned (if deriv=2)
  if(missing(deriv)){
    if(detection.function=="hazard") out <- 1 - exp(-(d/theta[1])^(-theta[2]))
    if(detection.function=="halfnormal") out <- exp(-d^2/(2*theta^2))
  }
  if(!missing(deriv)){
    if((deriv==1)&(detection.function=="halfnormal")) out <- d^2/theta[1]^3*detection.prob(d,"halfnormal",theta)
    if((deriv==1)&(detection.function=="hazard"))
      out <-  (detection.prob(d,"hazard",theta)-1) * (-1) * d / theta[1]^2 * theta[2] * (d/theta[1])^(-theta[2]-1)
    if((deriv==2)&(detection.function=="hazard"))
      out <-  (detection.prob(d,"hazard",theta)-1) * (d/theta[1])^(-theta[2]) * log(d/theta[1])
  }
  out
}
