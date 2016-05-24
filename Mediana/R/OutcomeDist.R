######################################################################################################################

# Function: OutcomeDist.
# Argument: Outcome Distribution and Outcome Type
# Description: This function is used to create an object of class Outcome.
#' @export
OutcomeDist = function(outcome.dist, outcome.type = NULL) {

  # Error checks
  if (!is.character(outcome.dist)) stop("Outcome: outcome distribution must be character.")
  if (!is.null(outcome.type)) {
    if (!is.character(unlist(outcome.type))) stop("Outcome: outcome must be character.")
    if(!all((unlist(outcome.type) %in% c("event","standard")))==TRUE) stop("Outcome: outcome type must be event or standard")
  }

  outcome = list(outcome.dist = outcome.dist,
                 outcome.type = outcome.type)

  class(outcome) = "OutcomeDist"
  return(outcome)
  invisible(outcome)

}