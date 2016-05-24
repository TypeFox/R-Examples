# Function: UniformDist.
# Argument: List of parameters (number of observations, maximum value).
# Description: This function is used to generate uniform outcomes.

UniformDist = function(parameter) {

  # Error checks
  if (missing(parameter))
    stop("Data model: UniformDist distribution: List of parameters must be provided.")

  if (is.null(parameter[[2]]$max))
    stop("Data model: UniformDist distribution: Maximum value must be specified.")

  max.value = parameter[[2]]$max

  if (max.value <= 0)
    stop("Data model: UniformDist distribution: Maximum value must be positive.")

  # Determine the function call, either to generate distribution or to return description
  call = (parameter[[1]] == "description")

  # Generate random variables
  if (call == FALSE) {
    # Error checks
    n = parameter[[1]]
    if (n%%1 != 0)
      stop("Data model: UniformDist distribution: Number of observations must be an integer.")
    if (n <= 0)
      stop("Data model: UniformDist distribution: Number of observations must be positive.")

    result = stats::runif(n = n, max = max.value)
  } else {
    # Provide information about the distribution function
    if (call == TRUE) {
      # Labels of distributional parameters
      result = list(list(max = "max"),list("Uniform"))
    }
  }
  return(result)

}
#End of UniformDist