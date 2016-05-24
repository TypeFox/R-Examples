#' Logit function
#'
#' Computes logit transformation.
#'
#' @param p probability
#' @return \code{ logit(p)} returns [log(p/(1-p)]
#' @author Jeff Laake
#' @keywords utility
logit <- function(p){
  log(p/(1-p))
}
