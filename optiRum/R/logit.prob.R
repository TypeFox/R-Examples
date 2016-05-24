#' Convert a logit to probability
#'
#' Transform a logit response from a glm into probability
#'
#' @param logit The log(odds)
#' @return prob Probability
#' 
#' @keywords logit odds glm probability
#' @seealso \code{\link{logit.odd}}  \code{\link{odd.prob}} 
#' @family creditrisk
#' @export
#' 
#' @examples
#' logit.prob(0) # equals 0.5

logit.prob <- function(logit) {
    odd.prob(logit.odd(logit))
} 
