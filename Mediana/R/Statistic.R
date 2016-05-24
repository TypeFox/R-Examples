######################################################################################################################

# Function: Statistic.
# Argument: Statistic ID, Statistical method, Samples and Parameters.
# Description: This function is used to create an object of class Statistic.
#' @export
Statistic = function(id, method, samples, par = NULL) {

  # Error checks
  if (!is.character(id)) stop("Statistic: ID must be character.")
  if (!is.character(method)) stop("Statistic: statistical method must be character.")
  if (!is.list(samples)) stop("Statistic: samples must be wrapped in a list.")
  if (any(lapply(samples, is.character) == FALSE)) stop("Statistic: samples must be character.")
  if (!is.null(par) & !is.list(par)) stop("MultAdj: par must be wrapped in a list.")

  statistic = list(id = id, method = method, samples = samples, par = par)

  class(statistic) = "Statistic"
  return(statistic)
  invisible(statistic)
}