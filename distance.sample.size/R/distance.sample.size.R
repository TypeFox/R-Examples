#' Required study size in distance sampling.
#'
#' Calculates the study size needed to achieve a target coefficient of variation for the abundance estimator in conventional distance sampling.
#' @param cv.pct The required cv expressed as a percentage. For example, use cv=15 for a coefficient of variation of 15\%.
#' @param N Optional. The total abundance of the objects or animals of interest in the whole region of interest. In practice may not be known, in which either
#' a rough estimate can be used, or N can be set to infinity (the default) which is equivalent to assuming that the fraction of all animals observed is small.
#' Setting N to Inf results in an over-estimation (usually slight) of the required sample size.
#' @param overdispersion The factor by which the variance of the number of objects observed is inflated due to overdispersion.
#'       Burnham, Anderson and Laake (1985) suggest that a value of 2 may be fairly typical in practice.
#' @param detection.function The detection function. Only "halfnormal" and "hazard" (hazard rate) are supported at present.
#' @param theta The detection function parameters. A single value for halfnormal, or a vector of two values for hazard rate.
#' @param mean.detection.prob.value An optional value specifying the mean detection probability over the range of observation. If this is supplied,
#'             the first element of theta should be set to NA, and theta[1] will be calculated using mean.detection.prob.value and detection.function.
#' @param shape.hazard Can be used to specify theta according to 3 preset hazard rate models (the ones used in the simulation in Clark 2016).
#'                     If shape.hazard is supplied, detection.function should be "hazard", theta need not be supplied, and w need not be supplied as is set to 1
#'                     (results in detection probabilities of 0.1 to 0.15 at w). All three options have an average detection rate of 0.6.
#' @param w The maximum range of observation. Objects at distance greater than w from the observer are assumed to never be recorded.
#' @param stop Set to T to open a browser window (for debugging purposes)
#' @details
#' It may be impossible to achieve the target precision, even if the expected sample size is equal to its maximum possible value of N divided by the mean detection probability.
#' In this case, missing values are returned for the required sample size and coverage proportion, and a warning is issued.
#' @return A vector with named values giving: the required expected sample size, the required coverage rate (i.e. the proportion P of the region falling within
#'         distance w of an observer's path), the penalty due to unknown detection parameters when P<<1, and the penalty due to unknown detection parameters for the
#'         required value of P. The user can then use either the required coverage rate to determine how closely to space transect lines (or how many points to
#'         select in a point transect study)
#' @references
#' Buckland S, Anderson D, Burnham K, Laake J and Borchers D (2001). Introduction to Distance Sampling: Estimating Abundance of Biological Populations. Oxford: Oxford University Press.
#'
#' Burnham KP, Anderson DR, Laake JL (1985), "Efficiency and bias in strip and line transect sampling". The Journal of Wildlife Management, pp. 1012-1018.
#'
#' Clark, R. G. (2016), "Statistical efficiency in distance sampling," PLoS One, forthcoming, www.plosone.org
#' @examples
#' distance.sample.size(cv.pct=15,N=1000,detection.function="hazard",shape.hazard="narrow")

distance.sample.size <- function(cv.pct,N=Inf,overdispersion=2,detection.function,theta,mean.detection.prob.value,shape.hazard=c("verynarrow","narrow","wide"),w,stop=F){
  if(stop) browser()
  cv <- cv.pct / 100
  # assign values to theta if shape.hazard supplied
  if(missing(theta)&(detection.function=="halfnormal")) theta <- NA
  if(missing(theta)&(detection.function=="hazard")) theta <- c(NA,NA)
  if((is.na(theta[2]))&(!missing(shape.hazard))){
    if(missing(mean.detection.prob.value)) mean.detection.prob.value <- 0.6
    if(shape.hazard=="verynarrow") theta[2] <- 1.25
    if(shape.hazard=="narrow") theta[2] <- 2
    if(shape.hazard=="wide") theta[2] <- 3
    w <- 1
  }
  # assign values to theta[1] if missing
  if(is.na(theta[1])&(!missing(mean.detection.prob.value))&(detection.function=="hazard"))
    theta[1] <- calculate.theta.given.mean.detection.prob(detection.function=detection.function,theta2=theta[2],w=w,mean.detection.prob.value=mean.detection.prob.value)
  if(is.na(theta[1])&(!missing(mean.detection.prob.value))&(detection.function=="halfnormal"))
    theta[1] <- calculate.theta.given.mean.detection.prob(detection.function=detection.function,w=w,mean.detection.prob.value=mean.detection.prob.value)
  # now calculate penalty0 (the penalty due to unknown detection parameters when P=0)
  penalty0 <- DS.penalty(detection.function=detection.function,theta=theta,w=w,P=0)
  # calculate mean detection probability if not specified
  if(missing(mean.detection.prob.value)) mean.detection.prob.value <- calculate.mean.detection.prob(theta=theta,w=w,detection.function=detection.function)
  # Now calculate required sample size
  required.sample.size <- penalty0 / ( cv^2/overdispersion + 1/N )
  if(required.sample.size>(N*mean.detection.prob.value)){
    required.sample.size <- NA
    warning("Not possible to achieve the target coefficient of variation. Note that even if the region is fully covered and all objects are detected, the distance sampling estimator has positive variance due to estimation of detection parameters.")
  }
  required.P <- required.sample.size / N
  penalty <- DS.penalty(detection.function=detection.function,theta=theta,w=w,P=required.P)
  if(N==Inf) required.P <- NA
  c( required.sample.size=required.sample.size , required.P=required.P , "penalty when sampling fraction assumed to be small"=penalty0 ,
     "actual penalty"=penalty )
}
