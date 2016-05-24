#' Variance penalty due to unknown detection parameters.
#'
#' Calculates the variance penalty factor due to unknown detection parameters in conventional distance sampling.
#' @param detection.function The detection function. Only "halfnormal" and "hazard" (hazard rate) are supported at present.
#' @param theta The detection function parameters. A single value for halfnormal, or a vector of two values for hazard rate.
#' @param mean.detection.prob.value An optional value specifying the mean detection probability over the range of observation. If this is supplied,
#'             the first element of theta should be set to NA, and theta[1] will be calculated using mean.detection.prob.value and detection.function.
#' @param w The maximum range of observation. Objects at distance greater than w from the observer are assumed to never be recorded.
#' @param P The proportion of the region of interest that is within w of an observer's path. P=0 may be assumed if the region is large relative to the observed area.
#' @param stop Set to T to open a browser window (for debugging purposes)
#' @return A single numeric value giving the asymptotic factor by which the variance is inflated due to unknown detection parameters.
#' @references
#' Buckland S, Anderson D, Burnham K, Laake J and Borchers D (2001). Introduction to Distance Sampling: Estimating Abundance of Biological Populations. Oxford: Oxford University Press.
#'
#' Clark, R. G. (2016), "Statistical efficiency in distance sampling," PLoS One, forthcoming, www.plosone.org
#' @examples
#' DS.penalty(detection.function="hazard",theta=c(NA,2),mean.detection.prob.value=0.6,w=1)

DS.penalty <- function(detection.function=c("halfnormal","hazard"),theta,mean.detection.prob.value,w,P=0,stop=F){
  if(stop) browser()
  if(is.na(theta[1])&(!missing(mean.detection.prob.value))&(detection.function=="hazard"))
    theta[1] <- calculate.theta.given.mean.detection.prob(detection.function=detection.function,theta2=theta[2],w=w,mean.detection.prob.value=mean.detection.prob.value)
  if(is.na(theta[1])&(!missing(mean.detection.prob.value))&(detection.function=="halfnormal"))
    theta[1] <- calculate.theta.given.mean.detection.prob(detection.function=detection.function,w=w,mean.detection.prob.value=mean.detection.prob.value)
  gbar.value <- calculate.mean.detection.prob(theta=theta,w,detection.function)
  hbar1 <- integrate(f=detection.prob,lower=0,upper=w,theta=theta,deriv=1,detection.function=detection.function)$value/w
  hbar2 <- 0
  if(length(theta)==2) hbar2 <- integrate(f=detection.prob,lower=0,upper=w,theta=theta,deriv=2,detection.function=detection.function)$value/w
  prodfn <- function(d,theta,detection.function,deriv1,deriv2){
    detection.prob(d=d,theta=theta,detection.function=detection.function,deriv=deriv1) *
      detection.prob(d=d,theta=theta,detection.function=detection.function,deriv=deriv2) /
      detection.prob(d=d,theta=theta,detection.function=detection.function)
  }
  Delta <- matrix(data=0,nrow=2,ncol=2)
  Delta[1,1] <- integrate( f=prodfn , lower=0 , upper=w , detection.function=detection.function , deriv1=1 , deriv2=1 , theta=theta)$value / (w*gbar.value) - (hbar1/gbar.value)^2
  if(length(theta)==2){
    Delta[2,2] <- integrate( f=prodfn , lower=0 , upper=w , detection.function=detection.function , deriv1=2 , deriv2=2 , theta=theta)$value / (w*gbar.value) - (hbar2/gbar.value)^2
    Delta[1,2] <- Delta[2,1] <- integrate( f=prodfn , lower=0 , upper=w , detection.function=detection.function , deriv1=1 , deriv2=2 , theta=theta)$value / (w*gbar.value) - hbar1*hbar2/gbar.value^2
  }
  hbar <- matrix( c( hbar1 , hbar2 ) , nrow=2,ncol=1)
  penalty <- 1 + t(hbar) %*% MASS::ginv(Delta) %*% hbar / gbar.value^2 / (1-P*gbar.value)
  as.numeric(penalty)
}

