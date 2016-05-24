#' Generalized logistic g function
#'
#' A parametric version of the g function following Richards (1959): 
#' \deqn{g(v) = M + \frac{1}{B} \log\left(\frac{Tv^T}{1-v^T}\right)}
#' @param v vector of standardized scores from the continuous ordinal scale, 0<v<1
#' @param par vector of 3 elements: \code{M}, the offset, \code{B}, the slope of the curve, and \code{T}, the symmetry of the curve 
#' @keywords Richards, generalized logistic function.
#' @details The generalized logistic functions maps from (0,1) to \eqn{(-\infty,\infty)}. 
#' \code{B} is the slope of the curve,  \code{T} is the symmetry and \code{M} is the offset. 
#' @return A vector of length equal to the length of \code{v}, with values \eqn{g(v)}.
#' 
#' @references Richards, F. (1959). A flexible growth function for empirical use, 
#' \emph{Journal of Experimental Botany}, 10, 290-301.
#' @seealso \code{\link{dg_glf}}, \code{\link{g_glf_inv}}
#' @author Maurizio Manuguerra, Gillian Heller
#


g_glf <- function(v, par){
  return(par[1] + log(par[3]*v^par[3]/(1-v^par[3]))/par[2])
}


#' Derivative of generalized logistic g function
#'
#' Derivative of the generalized logistic function as in Richards (1959): 
#' \deqn{g'(v) = \frac{T}{B}  \frac{1}{v(1-v^{T})}}
#' @param v vector of standardized scores from the continuous ordinal scale, 0<\code{v}<1.
#' @param par vector of 2 elements: \code{B}, the slope of the curve, and \code{T}, the symmetry of the curve.
#' @keywords Richards, derivative, generalized logistic function.
#' @return A vector of length equal to the length of \code{v}, with values \eqn{g'(v)}.
#' @seealso \code{\link{g_glf}},\code{\link{g_glf_inv}}
#' @author Maurizio Manuguerra, Gillian Heller
#'  @references Richards, F. (1959). A flexible growth function for empirical use, 
#' \emph{Journal of Experimental Botany}, 10, 290-301.
#

dg_glf <- function(v, par){
  #par = c(B, T)
  return(par[2]/par[1]/v/(1-v^par[2]))
}


#' Inverse of generalized logistic g function
#'
#' Inverse of a parametric version of the g function 
#' following Richards (1959): 
#' \deqn{g^{-1}(W) = \left( \frac{e^{B(W-M)}}{T+e^{B(W-M)}}  \right)^{\frac{1}{T}}} 
#' @param W vector of scores on the latent scale \eqn{(-\infty,\infty)}
#' @param par vector of 3 elements: \code{M}, the offset of the curve, \code{B}, the slope of the curve, and \code{T}, the symmetry 
#' of the curve 
#' @return A vector of length equal to the length of \code{W}, with values \eqn{g^{-1}(W)}
#' @keywords Richards, generalized logistic function.
#' @seealso \code{\link{g_glf}},\code{\link{dg_glf}}
#' @author Maurizio Manuguerra, Gillian Heller
#'  @references Richards, F. (1959). A flexible growth function for empirical use, 
#' \emph{Journal of Experimental Botany}, 10, 290-301.

g_glf_inv <- function(W, par){
  exp.part <- exp(par[2]*(W-par[1]))
  return((exp.part/(par[3]+exp.part))^(1/par[3]))
}

