######################################################################################################################

# Function: Test.
# Argument: Test ID, Statistical method, Samples and Parameters.
# Description: This function is used to create an object of class Test.
#' @export
Test = function(id, method, samples, par = NULL) {

  # Error checks
  if (!is.character(id)) stop("Test: ID must be character.")
  if (!is.character(method)) stop("Test: statistical method must be character.")
  if (!is.list(samples)) stop("Test: samples must be wrapped in a list.")
  if (all(lapply(samples, is.list) == FALSE) & any(lapply(samples, is.character) == FALSE)) stop("Test: samples must be character.")
  if (all(lapply(samples, is.list) == TRUE) & (!is.character(unlist(samples)))) stop("Test: samples must be character.")
  if (!is.null(par) & !is.list(par)) stop("Test: par must be wrapped in a list.")

  test = list(id = id, method = method, samples = samples, par = par)

  class(test) = "Test"
  return(test)
  invisible(test)
}
