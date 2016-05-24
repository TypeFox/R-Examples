################################################################################
# 
# OptimiseParametersThomas.R
# Version 1.0
# 30/01/2015
#
# Optimisation procedure of finding parameters that best fit observed data for
# the Thomas model
#
# Args:
#   Area: Vector of grain sizes for obsvered area of occupancies
#   Observed: Vector of observed area of occupancies
#   extent: Extent of study area (km2)
#   tolerance: tolerance for integration. Lower numbers allow for greater
#              accuracy but require longer processing times
#   model = "Thomas"
#
# Returns:
#   optim.pars: list of parameters estimated from optimisation procedure
#
################################################################################

OptimiseParametersThomas <- function(area, 
                                     observed, 
                                     extent, 
                                     tolerance = 1e-6,
                                     model = "Thomas",
                                     starting.params = NULL) {
  # Retrive residual function, downscaling function and starting parameters
  # for model of choice
  resid.fun <- getFunction(paste("Resid", model, sep = ""))
  pred.fun <- getFunction(paste("Predict", model, sep = ""))  
  if(is.null(starting.params)) {
    starting.pars <- get(paste("Params", model, sep = ""))
  }
  if(!is.null(starting.params)) {
    starting.pars <- starting.params
  }
  
  # Optimisation procedure
  optimisation <- minpack.lm::nls.lm(par = starting.pars,
                                     fn = resid.fun,
                                     tolerance = tolerance,
                                     area = area,
                                     extent = extent,
                                     observed = log(observed),
                                     control = minpack.lm::nls.lm.control(
                                       maxiter = 1000))
  optim.pars <- as.list(coef(optimisation))
  return(optim.pars)
}
