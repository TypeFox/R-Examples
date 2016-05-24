#-----------------------------------------------------------------------------#
# Compute the derivative(s) of a weight 
# 
# Takes the derivative of the \code{\link{wght_calc}} function with respect to each 
# parameter in \code{params}.
# 
# @param integrand the function to passed to the argument 'f' of \code{\link{integrate}},
# which is part of \code{\link{wght_calc}}.
# @param fixed.effects vector of fixed effect parameters
# @param random.effects vector random effect parameters
# @param allocation the allocation ratio for which to compute the weights
# @param hide.errors print \code{grad} error messages. Defaults to TRUE.
# @param ... additional arguments passed to integrand
# @return vector of derivatives with respect to element of params
# @export
# 
#-----------------------------------------------------------------------------#

wght_deriv_calc <- function(integrand,
                            fixed.effects,
                            random.effects = NULL,
                            allocation,
                            hide.errors = TRUE,
                            ...)
{  
  ## Necessary pieces ##
  integrand <- match.fun(integrand)
  dots <- list(...)
  
  ## Integrand and  arguments ##
  int.args <- append(get_args(integrand, dots),
                     get_args(integrate, dots))
  
  args <- append(append(int.args, get_args(numDeriv::grad, dots)),
                 list(func           = wght_calc, 
                      integrand      = integrand, 
                      allocation     = allocation,
                      fixed.effects  = fixed.effects,
                      random.effects = random.effects))
  
  params <- c(fixed.effects, random.effects)
  
  ## Compute Derivatives ##
  dervs <- sapply(1:length(params), function(i){
    args$x <- params[i]; args$pos <- i
    f <- try(do.call(numDeriv::grad, args = args), silent = hide.errors)
    return(ifelse(is(f, 'try-error'), NA, f))
  })
  return(dervs)
}
