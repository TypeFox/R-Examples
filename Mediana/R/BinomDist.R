######################################################################################################################

# Function: BinomDist .
# Argument: List of parameters (number of observations, proportion/probability of success).
# Description: This function is used to generate binomial outcomes (0/1).

BinomDist = function(parameter) {

  # Error checks
  if (missing(parameter))
    stop("Data model: BinomDist distribution: List of parameters must be provided.")


  if (is.null(parameter[[2]]$prop))
    stop("Data model: BinomDist distribution: Proportion must be specified.")


  prop = parameter[[2]]$prop


  if (prop < 0 | prop > 1)
    stop("Data model: BinomDist distribution: Proportion must be between 0 and 1.")


  # Determine the function call, either to generate distribution or to return description
  call = (parameter[[1]] == "description")

  # Generate random variables
  if (call == FALSE) {
    # Error checks
    n = parameter[[1]]
    if (n%%1 != 0)
      stop("Data model: BinomDist distribution: Number of observations must be an integer.")
    if (n <= 0)
      stop("Data model: BinomDist distribution: Number of observations must be positive.")


    result = stats::rbinom(n = n, size = 1, prob = prop)
  } else {
    # Provide information about the distribution function
    if (call == TRUE) {
      # Labels of distributional parameters
      result = list(list(prop = "prop"),list("Binomial"))
    }
  }
  return(result)

}
#End of BinomDist