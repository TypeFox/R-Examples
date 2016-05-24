#' Daily record of animal's movement (from 2012 to 2013).
#'
#' One dataset containing the number of animals that were moved from one node to 
#' another.
#'
#' @format A data frame with 78 rows and 5 variables:
#' \itemize{
#'   \item Finalidade: Type of destiny premises
#'   \item Dia: The day when the movement occurs
#'   \item originID: The ID of the origin premises
#'   \item destinyID: The ID of the destiny premises
#'   \item num.animais: The number of animals traded
#' }
#' @source ADAGRO
"networkSample"

#' Information about animal premises (from 2012 to 2013).
#'
#' A dataset containing animal premises' identification and census.
#'
#' @format A data frame with 507 rows and 2 variables:
#' \itemize{
#'   \item nodes.ID: The ID of the premises
#'   \item pop: premises's population size
#' }
"nodesCensus"