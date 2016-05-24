######################################################################################################################

# Function: MVExpoDist.
# Argument: List of parameters (number of observations, list(list(rate), correlation matrix).
# Description: This function is used to generate correlated exponential outcomes.

MVExpoDist = function(parameter) {

  # Error checks
  if (missing(parameter))
    stop("Data model: MVExpoDist distribution: List of parameters must be provided.")

  if (is.null(parameter[[2]]$par))
    stop("Data model: MVExpoDist distribution: Parameter list (rate) must be specified.")

  if (is.null(parameter[[2]]$corr))
    stop("Data model: MVExpoDist distribution: Correlation matrix must be specified.")

  par = parameter[[2]]$par
  corr = parameter[[2]]$corr


  # Number of endpoints
  m = length(par)


  if (ncol(corr) != m)
    stop("Data model: MVExpoDist distribution: The size of the hazard rate vector is different to the dimension of the correlation matrix.")
  if (sum(dim(corr) == c(m, m)) != 2)
    stop("Data model: MVExpoDist distribution: Correlation matrix is not correctly defined.")
  if (det(corr) <= 0)
    stop("Data model: MVExpoDist distribution: Correlation matrix must be positive definite.")
  if (any(corr < -1 | corr > 1))
    stop("Data model: MVExpoDist distribution: Correlation values must be comprised between -1 and 1.")


  # Determine the function call, either to generate distribution or to return description
  call = (parameter[[1]] == "description")


  # Generate random variables
  if (call == FALSE) {
    # Error checks
    n = parameter[[1]]
    if (n%%1 != 0)
      stop("Data model: MVExpoDist distribution: Number of observations must be an integer.")
    if (n <= 0)
      stop("Data model: MVExpoDist distribution: Number of observations must be positive.")

    # Generate multivariate normal variables
    multnorm = mvtnorm::rmvnorm(n = n, mean = rep(0, m), sigma = corr)
    # Store resulting multivariate variables
    mvmixed = matrix(0, n, m)
    # Convert selected components to a uniform distribution and then to exponential distribution
    for (i in 1:m) {
      uniform = stats::pnorm(multnorm[, i])

      if (is.null(par[[i]]$rate))
        stop("Data model: MVExpoDist distribution: Hazard rate parameter in the exponential distribution must be specified.")

      # Hazard rate
      hazard = as.numeric(par[[i]]$rate)
      if (hazard <= 0)
        stop("Data model: MVExpoDist distribution: Hazard rate parameter in the exponential distribution must be positive.")
      mvmixed[, i] = -log(uniform)/hazard
    }
    result = mvmixed

  } else {
    # Provide information about the distribution function
    if (call == TRUE) {
      # Labels of distributional parameters
      par.labels = list()
      for (i in 1:m) {
        par.labels[[i]] = list(rate = "rate")
      }
      result = list(list(par = par.labels, corr = "corr"),list("Multivariate Exponential"))
    }
  }
  return(result)
}
# End of MVExpoDist