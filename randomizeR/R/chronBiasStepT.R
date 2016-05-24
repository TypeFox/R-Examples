#' @include chronBias.R
NULL

###############################################
# --------------------------------------------#
# Class chronBiasStepT                        #
# --------------------------------------------#
###############################################

# Validity check function for objects of the chronBias class
# 
# @inheritParams overview 
#
# @return Returns a \code{TRUE}, if the settings of the object are valid.
validateChronBiasStepT <- function(object) {
  errors <- character()
  
  lengthSaltus<- length(object@saltus)
  if (lengthSaltus != 1) {
    msg <- paste("saltus is length ", lengthSaltus, ". Should be 1.", sep = "")
    errors <- c(errors, msg)
  }
  
  if(length(errors) == 0) TRUE else errors
}


# --------------------------------------------
# Class definition for chronBiasStepT
# --------------------------------------------

# Randomization paramters generic
setClass("chronBiasStepT",
         slots = c(saltus = "numeric"),
         contains = "chronBias",
         validity = validateChronBiasStepT)

# --------------------------------------------
# Constructor function for chronBiasStepT
# --------------------------------------------

# Generate chronBias object
# 
# The function generates an object of the S4 class \code{chronBias}. The object
# contains the information of possible chronological bias in a clinical trial.
# 
# @details
# The generated object contains full information
# \itemize{
# \item of the time trend function.
# \item of the strength of the time trend.
# \item of the saltus of the time trend.
# \item whether one test decision should be simulated
#   or the p.value should be calculated exact.
# \item of the alpha level of the two-sided test or of the 
#  quantiles of the corresponding distribution function used to determine 
#  an exact type-I-error probability of a given randomization sequence.
# }
# 
# @inheritParams chronBias
# @param type character string, should be "\code{stepT}":
# \describe{
#   \item{\code{Using "stepT" the following step function is used:}}{\deqn{f(i) 
#   = 1_{i \geq c >1} \theta}{1_{i \ge c >1} theta }}
# }
# 
# 
#
chronBiasStepT <- function(type, theta, method, saltus, alpha = 0.05) {
  new("chronBiasStepT", type = type, theta = theta, method = method, saltus = saltus, alpha = alpha)
}
