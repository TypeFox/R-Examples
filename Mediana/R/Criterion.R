######################################################################################################################

# Function: Criterion.
# Argument: Criterion ID, method, tests, statistics, parameters, labels.
# Description: This function is used to create an object of class Criterion
#' @export
Criterion = function(id, method, tests = NULL, statistics = NULL, par = NULL, labels) {

  # Error checks
  if (!is.character(id)) stop("Criterion: ID must be character.")
  if (!is.character(method)) stop("Criterion: method must be character.")
  if (!is.null(tests) & !is.list(tests)) stop("Criterion: tests must be wrapped in a list.")
  if (any(lapply(tests, is.character) == FALSE)) stop("Criterion: tests must be character.")
  if (!is.null(statistics) & !is.list(statistics)) stop("Criterion: statistics must be wrapped in a list.")
  if (any(lapply(statistics, is.character) == FALSE)) stop("Criterion: statistics must be character.")
  if (is.null(tests) & is.null(statistics )) stop("Criterion: tests and/or statistics must be provided")

  criterion = list(id = id ,
                   method = method ,
                   tests = tests ,
                   statistics = statistics ,
                   par = par ,
                   labels = labels)

  class(criterion) = "Criterion"
  return(criterion)
  invisible(criterion)
}