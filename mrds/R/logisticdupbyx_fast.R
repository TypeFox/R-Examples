#' Logistic for duplicates as a function of covariates (fast)
#'
#' As \code{\link{logisticdupbyx}}, but faster when distance is a covariate (but no interactions with distance occur.
#'
#' @inheritParams logisticdupbyx
#' @param beta_distance parameter for distance
#' @param x1 linear predictor for 1, without distance
#' @param x2 linear predictor for 2, without distance
#' @author David L Miller
logisticdupbyx_fast <- function(distance, x1, x2, models, beta, point, beta_distance){

  # function to calculate p/(1+p)
  ologit <- function(p) p/(1+p)

  # first part of the function
  gx1 <- ologit(exp(x1 + distance*beta_distance))

  # calculate second and return
  if(!point){
    return(gx1 * ologit(exp(x2 + distance*beta_distance)))
  }else{
    return(gx1 * ologit(exp(x2 + distance*beta_distance))*2*distance)
  }
}
