#' Convert an odd into a logit
#'
#' Transforming odds into logits (the response from binomial glms)
#'
#' @param odds Odds
#' @return logit Log(odds)
#'
#' @keywords logit odds glm probability
#' @seealso \code{\link{logit.odd}}  \code{\link{logit.prob}} 
#' @family creditrisk
#' @export
#' 
#' @examples
#' odd.logit(1) # equals 0

odd.logit <- function(odds) {
    log(odds)
} 
