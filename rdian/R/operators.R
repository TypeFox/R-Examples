#'@title AND and OR operators for Guardian filters and queries.
#'
#'@description \code{guardian_and} and \code{guardian_or} provide
#'(respectively) the AND and OR logical operators. If you pass them your
#'query terms, it passes them back either separated (so that the Guardian API
#'knows to consider a match to any \emph{one} term a match) or grouped (so that
#'the API only matches if every term appears).
#'
#'@param ... a vector of terms (or several vector of terms)
#'
#'@return a single string containing the terms, separated by the AND (,) or
#'OR (|) separators used by the Guardian API.
#'
#'@examples
#'# Simple AND example
#'guardian_and("sausage", "mash")
#'
#'# With ORs
#'guardian_or("sausage", "mash")
#'@aliases guardian_operators guardian_and guardian_or
#'@export
#'@rdname guardian_operators
guardian_and <- function(...){
  paste(..., sep = ",")
}

#'@rdname guardian_operators
#'@export
guardian_or <- function(...){
  paste(..., sep = "|")
} 
