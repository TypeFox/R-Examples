#' Function for internal use only.
#'
#' Calculates the value of the first element of the detection parameters theta based on a supplied value of the average detection rate.
#' Primarily for internal use.
#' @param detection.function The detection function. Only "halfnormal" and "hazard" (hazard rate) are supported at present.
#' @param theta2 The second detection parameter for the hazard rate model (the shape parameter). Not required for halfnormal.
#' @param w The maximum range of observation. Objects at distance greater than w from the observer are assumed to never be recorded.
#' @param mean.detection.prob.value The mean detection probability over the range of observation.
#' @param stop Set to T to open a browser window (for debugging purposes)
#' @return The value of the first parameter of theta.
#' @references
#' Buckland S, Anderson D, Burnham K, Laake J and Borchers D (2001). Introduction to Distance Sampling: Estimating Abundance of Biological Populations. Oxford: Oxford University Press.
#'
#' Clark, R. G. (2016), "Statistical efficiency in distance sampling," PLoS One, forthcoming, www.plosone.org
#' @examples
#' calculate.theta.given.mean.detection.prob(detection.function="hazard",theta2=2,w=1,
#' mean.detection.prob.value=0.6) # should equal 0.448
#' calculate.theta.given.mean.detection.prob(detection.function="halfnormal",w=1,
#' mean.detection.prob.value=0.6) # should equal 0.502

calculate.theta.given.mean.detection.prob <- function(detection.function,theta2,w,mean.detection.prob.value,stop=F){
  if(stop) browser()
  #integrand.fn <- function(d,theta)   detection.prob(d,detection.function=detection.function,theta=theta)
  mean.detection.prob.minus.target.value <- function(theta1,theta2=NULL,w,mean.detection.prob.value,detection.function=detection.function){
    calculate.mean.detection.prob(theta=c(theta1,theta2),w=w,detection.function=detection.function) - mean.detection.prob.value
  }
  if(detection.function=="halfnormal"){
    theta1 <- uniroot(f=mean.detection.prob.minus.target.value,interval=w*c(0.001,1000),w=w,
                      detection.function=detection.function,
                      mean.detection.prob.value=mean.detection.prob.value)$root
  }
  if(detection.function=="hazard"){
    # calculate search range such that g(w) is between 0.001 and 0.999
    theta1.range <- w * (-log(1-c(0.001,0.999)))^(1/theta2)
    theta1 <- uniroot(f=mean.detection.prob.minus.target.value,interval=theta1.range,
                      detection.function=detection.function,
                        w=w,mean.detection.prob.value=mean.detection.prob.value,theta2=theta2)$root

  }
  theta1
}
