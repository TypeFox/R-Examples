#' Data set on attachment dependency (Cook, 2000)

#' The classic Cook (2000) dataset consists of measurements on
#' security of attachment within families. Only the variable
#' measuring attachment dependency in family relationships is
#' included in this dataset. Four roles are present (i.e. two parents and two children):
#' mothers "m", fathers "f", the older
#' child as "c", and the younger child "y".
#'
#' The variables are as follows:
#'
#' \itemize{

#'   \item family.id An indicator for the family.
#'   \item actor.id An indicator for the perceiver, either "m", "f", "c", or "y". 
#'   \item partner.id An indicator for the target, either "m", "f", "c", or "y". 
#'   \item dep1 first measurement of attachment dependency. 
#'   \item dep2 second measurement of attachment dependency.
#' }
#'
#' @references Cook, W. L. (2000). Understanding attachment security in family context. \emph{Journal of Personality and Social Psychology, 78}, 285-294. doi:10.1037/0022-3514.78.2.285
#' @docType data
#' @keywords datasets
#' @format A data frame with 2496 rows and 4 variables (208 families with 4 members each, round-robin design)
#' @name two.indicators
#' @examples
#' data(two.indicators)
#' head(two.indicators)
NULL
