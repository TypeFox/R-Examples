#' Convert an odds to probability
#'
#' Transform odds into a probability
#'
#' @param odds Odds
#' @return prob Probability
#'
#' @keywords logit odds glm probability
#' @seealso \code{\link{logit.odd}}  \code{\link{logit.prob}} 
#' @family creditrisk
#' @export
#' 
#' @examples
#' odd.prob(1) # equals 0.5

odd.prob <- function(odds) {
    odds/(1 + odds)
} 
