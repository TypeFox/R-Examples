#' Data set on attachment anxiety in four person families (Cook, 2000)
#' 
#' The classic Cook (2000) dataset consists of measurements on security of attachment within families. 
#' Only the variable measuring fear of rejection in family relationships is included in this dataset. 
#' Cook (2000) deduced this variable from the anxiety dimension of the adult attachment scale (Collins & Read, 1990).
#' The original data consisted of four person families (i.e. two parents and two children): the mother is labeled as "m",  the father as "f", the oldest child as "c" and the youngest child as "y".
#'
#' The variables are as follows:
#' 
#' \itemize{
#'   \item family.id An indicator for the family.
#'   \item actor.id An indicator for the rater in the dyad, either "m", "f", "c", or "y"
#'   \item partner.id An indicator for the person being rated in the dyad, either "m", "f", "c", or "y"
#'   \item anx: The obtained score on the attachment anxiety scale.
#' }
#' 
#' @source This dataset was retrieved from \url{http://davidakenny.net/kkc/c9/c9.htm} in wide format and converted to an R dataset in long format.
#' @references Cook, W. L. (2000). Understanding attachment security in family context. \emph{Journal of Personality and Social Psychology, 78}, 285-294. doi:10.1037/0022-3514.78.2.285
#' @references Collins, N. L., & Read, S. J. (1990, April). Adult attachment, working models, and relationship quality in dating couples. Journal of Personality and Social Psychology, 58, 644-663.
#' @docType data
#' @keywords datasets
#' @format A data frame with 2497 rows and 4 variables (208 families with 4 members each, round-robin design)
#' @name four.person
#' @examples 
#' data(four.person)
#' head(four.person)
NULL
