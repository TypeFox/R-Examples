#' Data set on attachment anxiety in 3-person families (based on Cook, 2000)
#' 
#' The classic Cook (2000) dataset consists of measurements on security of attachment within families. 
#' Only the variable measuring fear of rejection in family relationships is included in this dataset. 
#' Cook (2000) deduced this variable from the anxiety dimension of the adult attachment scale (Collins & Read, 1990).
#' The orignial data consisted of four person families (i.e. two parents and two children), but in order to obtain  a three person family the oldest sibling is systematically deleted in accordance with Kenny, Kashy & Cook (2006).  The following labels are used for the three roles: the mother is labelled "m", the father "f" and the youngest child "y". The data are presented in the long format.
#' Three roles are present: Mothers "m", fathers "f", and the younger child "y".
#'
#' The variables are as follows:
#' 
#' \itemize{
#'   \item family.id An indicator for the family.
#'   \item actor.id An indicator of the rater in the dyad, either "m", "f", or "y"
#'   \item partner.id An indicator for the person being rated in the dyad, either "m", "f", or "y"
#'   \item anx1: First indicator of relationship specific anxiety (i.e. average of first halve of the scale).
#'   \item anx2: Second indicator of relationship specific anxiety (i.e. average of second halve of the scale).
#' }
#' @source This dataset was retrieved from \url{http://davidakenny.net/kkc/c9/c9.htm} in wide format and converted to an R dataset in long format.
#' @references Cook, W. L. (2000). Understanding attachment security in family context. \emph{Journal of Personality and Social Psychology, 78}, 285-294. doi:10.1037/0022-3514.78.2.285
#' @references Collins, N. L., & Read, S. J. (1990). Adult attachment, working models, and relationship quality in dating couples. \emph{Journal of Personality and Social Psychology}, 58, 644-663.
#' @docType data
#' @keywords datasets
#' @format A data frame with 1249 rows and 5 variables (208 families with 3 members each, round-robin design)
#' @name three.person
#' @examples 
#' data(three.person)
#' head(three.person)
NULL


