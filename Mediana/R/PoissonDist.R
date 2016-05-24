######################################################################################################################

# Function: PoissonDist .
# Argument: List of parameters (number of observations, mean).
# Description: This function is used to generate Poisson outcomes.

PoissonDist = function(parameter) {

  # Error checks
  if (missing(parameter))
    stop("Data model: PoissonDist distribution: List of parameters must be provided.")

  if (is.null(parameter[[2]]$lambda))
    stop("Data model: PoissonDist distribution: Lambda must be specified.")

  lambda = parameter[[2]]$lambda

  # Parameters check
  if (lambda <= 0)
    stop("Data model: PoissonDist distribution: Lambda must be non-negative.")

  # Determine the function call, either to generate distribution or to return description
  call = (parameter[[1]] == "description")

  # Generate random variables
  if (call == FALSE) {
    n = parameter[[1]]
    if (n%%1 != 0)
      stop("Data model: PoissonDist distribution: Number of observations must be an integer.")
    if (n <= 0)
      stop("Data model: PoissonDist distribution: Number of observations must be positive.")

    result = stats::rpois(n = n, lambda = lambda)
  } else {
    # Provide information about the distribution function
    if (call == TRUE) {
      # Labels of distributional parameters
      result = list(list(lambda = "lambda"),list("Poisson"))
    }
  }
  return(result)

}
# End of PoissonDist