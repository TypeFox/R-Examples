#' Logistic as a function of distance
#'
#' Treats logistic as a function of distance; for a given distance it computes
#' function at all covariate values in data.
#'
#' @param x covariate data
#' @param distance single distance value
#' @param models model list
#' @param beta logistic parameters
#'
#' @return vector of probabilities
#' @author Jeff Laake
logisticbyz <- function (x, distance, models, beta){

  #  Functions used: g0, setcov

  x$distance <- rep(distance, length(x$distance))
  zlist <- setcov(x, models$g0model)
  return(g0(beta,zlist))
}
