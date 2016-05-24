######################################################################################################################

# Function: SampleSize.
# Argument: A list or vector of numeric.
# Description: This function is used to create an object of class SampleSize.
#' @export
SampleSize = function(sample.size) {

  # Error checks
  if (any(!is.numeric(unlist(sample.size)))) stop("SampleSize: sample size must be numeric.")
  if (any(unlist(sample.size) %% 1 !=0)) stop("SampleSize: sample size must be integer.")
  if (any(unlist(sample.size) <=0)) stop("SampleSize: sample size must be strictly positive.")

  class(sample.size) = "SampleSize"
  return(unlist(sample.size))
  invisible(sample.size)

}