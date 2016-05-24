######################################################################################################################

# Function: MVMixedDist .
# Argument: List of parameters (number of observations, list(list (distribution type), list(distribution parameters) correlation matrix)).
# Description: This function is used to generate correlated normal, binary and exponential outcomes.

MVMixedDist = function(parameter) {

  # Error checks
  if (missing(parameter))
    stop("Data model: MVMixedDist distribution: List of parameters must be provided.")

  if (is.null(parameter[[2]]$type))
    stop("Data model: MVMixedDist distribution: Distribution type must be specified.")

  if (is.null(parameter[[2]]$par))
    stop("Data model: MVMixedDist distribution: Parameters list must be specified.")

  if (is.null(parameter[[2]]$corr))
    stop("Data model: MVMixedDist distribution: Correlation matrix must be specified.")


  type = parameter[[2]]$type
  par = parameter[[2]]$par
  corr = parameter[[2]]$corr


  # Number of endpoints
  m = length(par)


  if (length(type) != m)
    stop("Data model: MVMixedDist distribution: Number of distribution type parameters must be equal to the number of endpoints.")
  for (i in 1:m) {
    if ((type[[i]] %in% c("NormalDist", "BinomDist", "ExpoDist")) == FALSE)
      stop("Data model: MVMixedDist distribution: MVMixedDist accepts only normal, binomial and exponential endpoints.")
  }
  if (ncol(corr) != m)
    stop("Data model: MVMixedDist distribution: The size of the outcome parameter is different to the dimension of the correlation matrix.")
  if (sum(dim(corr) == c(m, m)) != 2)
    stop("Data model: MVMixedDist distribution: Correlation matrix is not correctly defined.")
  if (det(corr) <= 0)
    stop("Data model: MVMixedDist distribution: Correlation matrix must be positive definite.")
  if (any(corr < -1 | corr > 1))
    stop("Data model: MVMixedDist distribution: Correlation values must be comprised between -1 and 1.")


  # Determine the function call, either to generate distribution or to return description
  call = (parameter[[1]] == "description")


  # Generate random variables
  if (call == FALSE) {
    # Error checks
    n = parameter[[1]]
    if (n%%1 != 0)
      stop("Data model: MVMixedDist distribution: Number of observations must be an integer.")
    if (n <= 0)
      stop("Data model: MVMixedDist distribution: Number of observations must be positive.")

    # Generate multivariate normal variables
    multnorm = mvtnorm::rmvnorm(n = n, mean = rep(0, m), sigma = corr)

    # Store resulting multivariate variables
    mvmixed = matrix(0, n, m)

    # Convert selected components to a uniform distribution and then to either binomial or exponential distribution
    for (i in 1:m) {
      if (type[[i]] == "NormalDist") {
        if (is.null(par[[i]]$mean))
          stop("Data model: MVMixedDist distribution: Mean in the normal distribution must be specified.")
        if (is.null(par[[i]]$sd))
          stop("Data model: MVMixedDist distribution: SD in the normal distribution must be specified.")
        mean = as.numeric(par[[i]]$mean)
        sd = as.numeric(par[[i]]$sd)
        if (sd <= 0)
          stop("Data model: MVMixedDist distribution: SD in the normal distribution must be positive.")
        mvmixed[, i] = mean + sd * multnorm[, i]
      } else if (type[[i]] == "BinomDist") {
        uniform = stats::pnorm(multnorm[, i])
        # Proportion
        if (is.null(par[[i]]$prop))
          stop("Data model: MVMixedDist distribution: Proportion in the binomial distribution must be specified.")
        prop = as.numeric(par[[i]]$prop)
        if (prop < 0 | prop > 1)
          stop("Data model: MVMixedDist distribution: Proportion in the binomial distribution must be between 0 and 1.")
        mvmixed[, i] = (uniform <= prop)
      } else if (type[[i]] == "ExpoDist") {
        uniform = stats::pnorm(multnorm[, i])
        # Hazard rate
        if (is.null(par[[i]]$rate))
          stop("Data model: MVMixedDist distribution: Hazard rate in the exponential distribution must be specified.")
        hazard = as.numeric(par[[i]])
        if (hazard <= 0)
          stop("Data model: MVMixedDist distribution: Hazard rate parameter in the exponential distribution must be positive.")
        mvmixed[, i] = -log(uniform)/hazard
      }
    }
    result = mvmixed

  } else {
    # Provide information about the distribution function
    if (call == TRUE) {
      # Labels of distributional parameters
      par.labels = list()
      outcome.name=""
      for (i in 1:m) {
        if (type[[i]] == "NormalDist")
        {
          par.labels[[i]] = list(mean = "mean", sd = "SD")
          outcome.name=paste0(outcome.name,", ","Normal")
        }
        if (type[[i]] == "BinomDist")
        {
          par.labels[[i]] = list(prop = "prop")
          outcome.name=paste0(outcome.name,", ","Binomial")
        }
        if (type[[i]] == "ExpoDist")
        {
          par.labels[[i]] = list(rate = "rate")
          outcome.name=paste0(outcome.name,", ","Exponential")
        }
      }
      result = list(list(type = "type", par = par.labels, corr = "corr"),list(paste0("Multivariate Mixed (", sub(", ","",outcome.name),")")))
    }
  }
  return(result)
}
# End of MVMixedDist