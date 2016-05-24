#' Logistic for duplicates as a function of covariates
#'
#' Treats logistic for duplicates as a function of covariate z; for a given z
#' it computes the function at with those covariate values at a range of
#' distances.
#'
#' @param distance vector of distance values
#' @param x1 covariate data for fct 1
#' @param x2 covariate data for fct 2
#' @param models model list
#' @param beta logistic parameters
#' @param point \code{TRUE} for point transect data
#'
#' @return vector of probabilities
#' @author Jeff Laake
logisticdupbyx <- function(distance, x1, x2, models, beta, point){

  # avoid using g0 which calls exp and matrix multiplication twice
  ologit <- function(p) p/(1+p)

  #  Functions used: g0, setcov
  xlist <- as.list(x1)
  xlist$distance <- distance
  xmat <- expand.grid(xlist)

  gx1 <- ologit(exp(setcov(xmat, models$g0model) %*% beta))

  xlist <- as.list(x2)
  xlist$distance <- distance
  xmat <- expand.grid(xlist)

  if(!point){
    return(gx1 * ologit(exp(setcov(xmat, models$g0model) %*% beta)))
  }else{
    return(gx1 * ologit(exp(setcov(xmat, models$g0model) %*% beta))*2*distance)
  }
}
