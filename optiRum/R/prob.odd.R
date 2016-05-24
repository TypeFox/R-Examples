#' Convert a probability into odds probability
#'
#' Transform probabilities into odds
#'
#' @param prob Probability
#' @return odds Odds
#'
#' @keywords logit odds glm probability
#' @seealso \code{\link{prob.logit}}  \code{\link{odd.logit}} 
#' @family creditrisk
#' @export
#' 
#' @examples
#' prob.odd(0.5) # equals 1

prob.odd <- function(prob) {
    prob/(1 - prob)
} 
