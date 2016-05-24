################################################################################
# 
# OptimiseParameters.R
# Version 1.0
# 28/01/2015
#
# Optimisation procedure of finding parameters that best fit observed data
#
# Args:
#   Area: Vector of grain sizes for obsvered area of occupancies
#   Observed: Vector of observed area of occupancies
#   Model: which downscaling model to use. Choice of:
#     Nachman   Nachman model
#     PL        Power Law model
#     Poisson   Poisson model
#     NB        Negative binomial model
#     GNB       Generalised negative binomial model
#     INB       Improved negative binomial model
#
# Returns:
#   optim.pars: list of parameters estimated from optimisation procedure
#
################################################################################

OptimiseParameters <- function(area, observed, model, starting.params = NULL) {
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
                                     control = minpack.lm::nls.lm.control(
                                       maxiter = 1000))
  optim.pars <- as.list(coef(optimisation))
  return(optim.pars)
}
