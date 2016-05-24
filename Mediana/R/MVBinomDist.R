######################################################################################################################

# Function: MVBinomDist .
# Argument: List of parameters (number of observations, list(list (prop), correlation matrix)).
# Description: This function is used to generate correlated multivariate binomial (0/1) outcomes.

MVBinomDist = function(parameter) {

  if (missing(parameter))
    stop("Data model: MVBinomDist distribution: List of parameters must be provided.")
  # Error checks
  if (is.null(parameter[[2]]$par))
    stop("Data model: MVBinomDist distribution: Parameter list (prop) must be specified.")
  if (is.null(parameter[[2]]$par))
    stop("Data model: MVBinomDist distribution: Correlation matrix must be specified.")

  par = parameter[[2]]$par
  corr = parameter[[2]]$corr

  # Number of endpoints
  m = length(par)

  if (ncol(corr) != m)
    stop("Data model: MVBinomDist distribution: The size of the proportion vector is different to the dimension of the correlation matrix.")
  if (sum(dim(corr) == c(m, m)) != 2)
    stop("Data model: MVBinomDist distribution: Correlation matrix is not correctly defined.")
  if (det(corr) <= 0)
    stop("Data model: MVBinomDist distribution: Correlation matrix must be positive definite.")
  if (any(corr < -1 | corr > 1))
    stop("Data model: MVBinomDist distribution: Correlation values must be comprised between -1 and 1.")


  # Determine the function call, either to generate distribution or to return description
  call = (parameter[[1]] == "description")


  # Generate random variables
  if (call == FALSE) {
    # Error checks
    n = parameter[[1]]
    if (n%%1 != 0)
      stop("Data model: MVBinomDist distribution: Number of observations must be an integer.")
    if (n <= 0)
      stop("Data model: MVBinomDist distribution: Number of observations must be positive.")

    # Generate multivariate normal variables
    multnorm = mvtnorm::rmvnorm(n = n, mean = rep(0, m), sigma = corr)
    # Store resulting multivariate variables
    mvbinom = matrix(0, n, m)
    # Convert selected components to a uniform distribution and then to binomial distribution
    for (i in 1:m) {
      uniform = stats::pnorm(multnorm[, i])
      # Proportion
      if (is.null(par[[i]]$prop))
        stop("Data model: MVBinomDist distribution: Proportion must be specified.")

      prop = as.numeric(par[[i]]$prop)
      if (prop < 0 | prop > 1)
        stop("Data model: MVBinomDist distribution: proportion in the binomial distribution must be between 0 and 1.")

      mvbinom[, i] = (uniform <= prop)
    }
    result = mvbinom

  } else {
    # Provide information about the distribution function
    if (call == TRUE) {
      # Labels of distributional parameters
      par.labels = list()
      for (i in 1:m) {
        par.labels[[i]] = list(prop = "prop")
      }
      result = list(list(par = par.labels, corr = "corr"),list("Multivariate Binomial"))
    }
  }
  return(result)
}
#End of MVBinomDist