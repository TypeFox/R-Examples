######################################################################################################################

# Function: Interim.
# Argument: Sample (sample), Criterion (criterion) and Fraction (fraction).
# Description: This function is used to create an object of class MultAdj.

Interim = function(sample = NULL, criterion = NULL, fraction = NULL) {

  if (is.null(sample))
    stop("Interim: a sample must be defined")
  if (is.null(criterion))
    stop("Interim: a criterion must be defined")
  if (is.null(fraction))
    stop("Interim: a fraction must be defined")

  interim = list(sample = sample, criterion = criterion ,fraction = fraction )

  class(interim) = "Interim"
  return(interim)
  invisible(interim)

}