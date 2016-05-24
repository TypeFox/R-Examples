# Function: NormalDist.
# Argument: List of parameters (number of observations, list(mean, standard deviation)).
# Description: This function is used to generate normal outcomes.
NormalDist = function(parameter) {

  # Error checks
  if (missing(parameter))
    stop("Data model: NormalDist distribution: List of parameters must be provided.")

  if (is.null(parameter[[2]]$mean))
    stop("Data model: NormalDist distribution: Mean must be specified.")

  if (is.null(parameter[[2]]$sd))
    stop("Data model: NormalDist distribution: SD must be specified.")

  mean = parameter[[2]]$mean
  sd = parameter[[2]]$sd

  if (sd <= 0)
    stop("Data model: NormalDist distribution: Standard deviations in the normal distribution must be positive.")

  # Determine the function call, either to generate distribution or to return description
  call = (parameter[[1]] == "description")

  # Generate random variables
  if (call == FALSE) {
    # Error checks
    n = parameter[[1]]
    if (n%%1 != 0)
      stop("Data model: NormalDist distribution: Number of observations must be an integer.")
    if (n <= 0)
      stop("Data model: NormalDist distribution: Number of observations must be positive.")

    result = stats::rnorm(n = n, mean = mean, sd = sd)

  } else {
    # Provide information about the distribution function
    if (call == TRUE) {
      # Labels of distributional parameters
      result = list(list(mean = "mean", sd = "SD"),list("Normal"))
    }
  }
  return(result)
}
#End of NormalDist