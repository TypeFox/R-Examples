######################################################################################################################

# Function: Sample.
# Argument: Sample ID, Outcome Parameters, Sample Size.
# Description: This function is used to create an object of class Sample.
#' @export
Sample = function(id, outcome.par,  sample.size = NULL) {

  # Error checks
  if (!is.character(unlist(id))) stop("Sample: sample ID must be character.")
  if (!is.list(outcome.par)) stop("Sample: outcome parameters must be provided in a list.")

  if (!is.null(sample.size)){
    # Error checks
    if (any(!is.numeric(unlist(sample.size)))) stop("Sample: sample size must be numeric.")
    if (any(unlist(sample.size) %% 1 !=0)) stop("Sample: sample size must be integer.")
    if (any(unlist(sample.size) <=0)) stop("Sample: sample size must be strictly positive.")

  }

  sample = list(id = id,
                outcome.par = outcome.par,
                sample.size = sample.size)

  class(sample) = "Sample"
  return(sample)
  invisible(sample)
}
