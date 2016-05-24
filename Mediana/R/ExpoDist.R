######################################################################################################################

# Function: ExpoDist.
# Argument: List of parameters (number of observations, rate).
# Description: This function is used to generate exponential outcomes.

ExpoDist = function(parameter) {

  # Error checks
  if (missing(parameter))
    stop("Data model: ExpoDist distribution: List of parameters must be provided.")

  if (is.null(parameter[[2]]$rate))
    stop("Data model: ExpoDist distribution: Rate parameter must be specified.")

  rate = parameter[[2]]$rate

  # Parameters check
  if (rate <= 0) stop("Data model: ExpoDist distribution: Rate parameter must be positive")

  # Determine the function call, either to generate distribution or to return description
  call = (parameter[[1]] == "description")

  # Generate random variables
  if (call == FALSE) {
    n = parameter[[1]]
    if (n%%1 != 0)
      stop("Data model: ExpoDist distribution: Number of observations must be an integer.")
    if (n <= 0)
      stop("Data model: ExpoDist distribution: Number of observations must be positive.")

    result = stats::rexp(n = n, rate = rate)
  } else {
    # Provide information about the distribution function
    if (call == TRUE) {
      # Labels of distributional parameters
      result = list(list(rate = "rate"),list("Exponential"))
    }
  }
  return(result)
}
# End of ExpoDist