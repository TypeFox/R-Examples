#' Convert a logit to odds
#'
#' Transform a logit response from a glm into odds
#'
#' @param logit The log(odds)
#' @return odds Odds
#' 
#' @keywords logit odds glm
#' @seealso \code{\link{logit.prob}} 
#' @family creditrisk
#' @export
#' 
#' @examples
#' logit.odd(0) # equals 1

logit.odd <- function(logit) {
    exp(logit)
} 
