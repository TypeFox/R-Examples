#' @name cv
#' @aliases cv
#' @author Lucas Venezian Povoa \email{lucasvenez@@gmail.com}
#' 
#' @title Coefficient of Variance
#' @description Calculates the coefficient of variance of a monthly precipitation.
#' @details A daily precipitation serie is transformed to a monthly serie.
#' @usage cv(object)
#' @param object is a daily or monthly precipitation serie.
#' @return the coefficient of variance
#' @examples 
#' ##
#' # Loading the montly precipitation serie.
#' data(monthly)
#' 
#' ##
#' # Calculating the Coefficient of Variance
#' cv(monthly)
#' @export
cv <- function(object) {
  object <- as.monthly(object)
  return(sd(object$precipitation)/mean(object$precipitation))
}