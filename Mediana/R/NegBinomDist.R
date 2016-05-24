######################################################################################################################

# Function: NegBinomDist .
# Argument: List of parameters (number of observations, list(dispersion, mean)).
# Description: This function is used to generate negative-binomial outcomes.

NegBinomDist = function(parameter) {

  # Error checks
  if (missing(parameter))
    stop("Data model: NegBinomDist distribution: List of parameters must be provided.")

  if (is.null(parameter[[2]]$dispersion))
    stop("Data model: NegBinomDist distribution: Dispersion (size) must be specified.")

  if (is.null(parameter[[2]]$mean))
    stop("Data model: NegBinomDist distribution: Mean (mu) must be specified.")

  dispersion = parameter[[2]]$dispersion
  mean = parameter[[2]]$mean

  # Parameters check
  if (dispersion <= 0) {
    stop("Data model: NegBinomDist distribution: Dispersion parameter must be positive.")
  } else if (mean <= 0) {
    stop("Data model: NegBinomDist distribution: Mean must be positive.")
  }

  # Determine the function call, either to generate distribution or to return description
  call = (parameter[[1]] == "description")

  # Generate random variables
  if (call == FALSE) {
    n = parameter[[1]]
    if (n%%1 != 0)
      stop("Data model: NegBinomDist distribution: Number of observations must be an integer.")
    if (n <= 0)
      stop("Data model: NegBinomDist distribution: Number of observations must be positive.")

    result = stats::rnbinom(n = n, size = dispersion, mu = mean)
  } else {
    # Provide information about the distribution function
    if (call == TRUE) {
      # Labels of distributional parameters
      result = list(list(dispersion = "dispersion", mean = "mean"),list("Negative binomial"))
    }
  }
  return(result)
}
#End of NegBinomDist