#' Produce a scaled score based on a logit
#'
#' This function takes a logit and scales based on an intercept and doubling of odds ratio
#'
#' @param logit Logit to be scaled
#' @param offset Midrange, default is 300
#' @param scale Value in which odds are double, default is 50
#' 
#' @keywords glm logit score
#' @seealso \code{\link{glm}}
#' @family creditrisk
#' @export
#' 
#' @examples
#' scaledScore(0) # 300
#' scaledScore(0,offset=600) # 600

scaledScore <- function(logit, offset = 300, scale = 50) {
    round(offset + logit * scale/log(2))
} 
