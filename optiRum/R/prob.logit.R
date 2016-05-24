#' Convert a probability into a logit
#'
#' Transforming probabilities into logits (the response from binomial glms)
#'
#' @param prob Probability
#' @return logit Log(odds)
#'
#' @keywords logit odds glm probability
#' @seealso \code{\link{prob.odd}}  \code{\link{odd.logit}} 
#' @family creditrisk
#' @export
#' 
#' @examples
#' prob.logit(0.5) # equals 0

prob.logit <- function(prob) {
    odd.logit(prob.odd(prob))
} 
