#'Extract coefficients from a fitted Stochastic Mortality Model
#'
#' @param object an object of class \code{"fitStMoMo"} with the fitted 
#' parameters of a stochastic mortality model.
#' @param ... other arguments.
#' 
#' @return A list of model parameters with components:
#'   
#'   \item{ax}{ Vector with the fitted values of the static age function 
#'   \eqn{\alpha_x}. If the model does not have a static age function or failed 
#'   to fit this is set to \code{NULL}.}
#'     
#'   \item{bx}{ Matrix with the values of the period age-modulating functions 
#'   \eqn{\beta_x^{(i)}, i=1, ..., N}. If the \eqn{i}-th age-modulating 
#'   function is non-parametric (e.g., as in the Lee-Carter model) 
#'   \code{bx[, i]} contains the estimated values. If the model does not have 
#'   any age-period terms (i.e. \eqn{N=0}) or failed to fit this is set to 
#'   \code{NULL}.}
#'   
#'   \item{kt}{ Matrix with the values of the fitted period indexes 
#'   \eqn{\kappa_t^{(i)}, i=1, ..., N}. \code{kt[i, ]} contains the estimated 
#'   values of the \eqn{i}-th period index. If the model does not have any 
#'   age-period terms (i.e. \eqn{N=0}) or failed to fit this is set to 
#'   \code{NULL}.}
#'   
#'   \item{b0x}{ Vector with the values of the cohort age-modulating function 
#'   \eqn{\beta_x^{(0)}}. If the age-modulating function is non-parametric 
#'   \code{b0x} contains the estimated values. If the model does not have a 
#'   cohort effect or failed to fit this is set to \code{NULL}.}
#'     
#'   \item{gc}{ Vector with the fitted cohort index \eqn{\gamma_{c}}.
#'   If the model does not have a cohort effect or failed to fit this is set 
#'   to \code{NULL}.}
#'   
#' @examples
#' APCfit <- fit(apc(), Dxt = EWMaleData$Dxt, Ext = EWMaleData$Ext, 
#'              ages = EWMaleData$ages, years = EWMaleData$years)
#' coef(APCfit) 
#' @export
coef.fitStMoMo <- function(object, ...) {
  list(ax = object$ax, bx = object$bx, kt = object$kt, b0x = object$b0x, 
       gc = object$gc)
}