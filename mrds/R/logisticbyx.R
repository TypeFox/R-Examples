#' Logistic as a function of covariates
#'
#' treats logistic as a function of covariates; for a given covariate
#' combination it computes function at with those covariate values at a
#' range of distances
#'
#' @param distance vector of distance values
#' @param x covariate data
#' @param models model list
#' @param beta logistic parameters
#' @param point \code{TRUE} if a point transect model
#'
#' @return vector of probabilities
#' @author Jeff Laake
logisticbyx <- function (distance, x, models, beta, point){

  # Functions used: g0, setcov

  xlist <- as.list(x)
  xlist$distance <- distance
  xmat <- expand.grid(xlist)

  if(!point){
    return(g0(beta, setcov(xmat, models$g0model)))
  }else{
    return(g0(beta, setcov(xmat, models$g0model))*2*distance)
  }
}
