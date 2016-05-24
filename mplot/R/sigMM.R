#' Standard deviation of hatQ_M and hatQ_tildeM
#' 
#' This function calculates \eqn{\hat(sigma)_{M,tildeM}}
#'  the standard deviation of
#' the difference between the two lack-of-fit measures
#' \eqn{\hat{Q}_M} and \eqn{\hat{Q}_{\tilde{M}}} as 
#' described in Jiang et. al. (2008).  When using the
#' adaptive fence procedure, this quantity no longer needs
#' to be calculated and simply returns a value of 1.
#' 
#' @param k.mod number of parameters in the estimated model
#' @param method method the model selection method to be used. Currently
#'   only \code{method = "ML"} is supported (perhaps in the future
#'   \code{method = "MVC"} will be implemented).
#' @param k.full number of parameters in the full model
#' @param adaptive logical. If \code{TRUE} the boundary of the fence is 
#'   given by cstar.  Otherwise, it the original (non-adaptive) fence
#'   is performed where the boundary is cstar*hat(sigma)_{M,tildeM}.
#' @noRd
sigMM = function(k.mod,method,k.full,adaptive){
  if(method=="ML" & adaptive==FALSE){
    return(sqrt((k.full-k.mod)/2))
  } else if (adaptive==TRUE){
    return(1)
  }
}