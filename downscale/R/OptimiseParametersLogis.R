################################################################################
# 
# OptimiseParametersLogis.R
# Version 1.0
# 05/02/2015
#
# Optimisation procedure of finding parameters that best fit observed data for
# the Logistic model
#
# Args:
#   Area: Vector of grain sizes for obsvered area of occupancies
#   Observed: Vector of observed area of occupancies
#   model = "Logis"
#
# Returns:
#   optim.pars: list of parameters estimated from optimisation procedure
#
################################################################################

OptimiseParametersLogis <- function(area, 
                                    observed,
                                    model = "Logis",
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
                                     area = area,
                                     observed = log(observed),
                                     lower = c("C" = 0, "z" = -Inf),
                                     upper = c("C" = Inf, "z" = Inf),
                                     control = minpack.lm::nls.lm.control(
                                       maxiter = 1000))
  optim.pars <- as.list(coef(optimisation))
  return(optim.pars)
}
