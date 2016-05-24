#' Threshold key function
#'
#' @param distance perpendicular distance vector
#' @param key.scale vector of scale values
#' @param key.shape vector of shape values
#'
#' @return vector of probabilities
keyfct.th1 <- function(distance, key.scale, key.shape){
  return( 0.5 - 0.5*erf(distance/key.scale - key.shape))
}


#' Threshold key function
#'
#' @param distance perpendicular distance vector
#' @param key.scale vector of scale values
#' @param key.shape vector of shape values
#'
#' @return vector of probabilities
keyfct.th2 <- function(distance, key.scale, key.shape){
  return(erf(exp(key.shape-distance/key.scale)))
}


# error function
erf <- function(x) 2*pnorm(x*sqrt(2)) - 1
