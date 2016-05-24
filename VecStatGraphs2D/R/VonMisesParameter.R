VonMisesParameter <- function (azimuths) 
{

  # Estimate the Von Mises concentration parameter for circular data
  # The method uses the approximation described in 
  # Fisher "Statistical Analysis of Circular Data"
  #
  # Args:
  #   azimuths: azimuths (degrees)
  #
  # Returns:
  #   Von Mises parameter

  n_elements = length(azimuths)
  mean_module = (MeanModule(azimuths))
  if (mean_module < 0.53) {
    parameter = (2 * mean_module) + (mean_module^3) + (5 * (mean_module^5)/6)
  }
  else{
	if ((mean_module >= 0.53) & (mean_module < 0.85)) {
		parameter = -0.4 + (1.39 * mean_module) + (0.43/(1 - mean_module))
	}
	else{
		# if (mean_module >= 0.85)
		parameter = 1/((mean_module^3) - (4 * mean_module^2) + 3 * mean_module)
	}
  }
  # correction for small samples 

  if (n_elements < 16) {
    if (parameter < 2) {
      parameter = max(parameter - (2 / (n_elements * parameter)), 0)
    }
    else  {
      parameter = (n_elements - 1)^3 * parameter / (n_elements^3 + n_elements)
    }
  }

  return(parameter)

}
