######################################################################################################################

# Function: MVNormalDist.
# Argument: List of parameters (number of observations, list(list (mean, SD), correlation matrix)).
# Description: This function is used to generate correlated multivariate normal outcomes.

MVNormalDist = function(parameter) {

  # Error checks
  if (missing(parameter))
    stop("Data model: MVNormalDist distribution: List of parameters must be provided.")

  if (is.null(parameter[[2]]$par))
    stop("Data model: MVNormalDist distribution: Parameter list (means and SDs) must be specified.")

  if (is.null(parameter[[2]]$corr))
    stop("Data model: MVNormalDist distribution: Correlation matrix must be specified.")

  par = parameter[[2]]$par
  corr = parameter[[2]]$corr

  # Number of endpoints
  m = length(par)

  if (ncol(corr) != m)
    stop("Data model: MVNormalDist distribution: The size of the mean vector is different to the dimension of the correlation matrix.")
  if (sum(dim(corr) == c(m, m)) != 2)
    stop("Data model: MVNormalDist distribution: Correlation matrix is not correctly defined.")
  if (det(corr) <= 0)
    stop("Data model: MVNormalDist distribution: Correlation matrix must be positive definite.")
  if (any(corr < -1 | corr > 1))
    stop("Data model: MVNormalDist distribution: Correlation values must be between -1 and 1.")


  # Determine the function call, either to generate distribution or to return description
  call = (parameter[[1]] == "description")


  # Generate random variables
  if (call == FALSE) {
    # Error checks
    n = parameter[[1]]
    if (n%%1 != 0)
      stop("Data model: MVNormalDist distribution: Number of observations must be an integer.")
    if (n <= 0)
      stop("Data model: MVNormalDist distribution: Number of observations must be positive.")


    # Generate multivariate normal variables
    multnorm = mvtnorm::rmvnorm(n = n, mean = rep(0, m), sigma = corr)
    # Store resulting multivariate variables
    mv = matrix(0, n, m)


    for (i in 1:m) {
      if (is.null(par[[i]]$mean))
        stop("Data model: MVNormalDist distribution: Mean must be specified.")
      if (is.null(par[[i]]$sd))
        stop("Data model: MVNormalDist distribution: SD must be specified.")

      mean = as.numeric(par[[i]]$mean)
      sd = as.numeric(par[[i]]$sd)

      if (sd <= 0)
        stop("Data model: MVNormalDist distribution: Standard deviations must be positive.")

      mv[, i] = mean + sd * multnorm[, i]
    }
    result = mv
  } else {
    # Provide information about the distribution function
    if (call == TRUE) {
      # Labels of distributional parameters
      par.labels = list()
      for (i in 1:m) {
        par.labels[[i]] = list(mean = "mean", sd = "SD")
      }
      result = list(list(par = par.labels, corr = "corr"),list("Multivariate Normal"))
    }
  }
  return(result)
}
#End of MVNormalDist