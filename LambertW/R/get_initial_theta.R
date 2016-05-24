#' @rdname theta-utils
#' 
#' @description
#' \code{get_initial_theta} provides initial estimates for \eqn{\alpha},
#'     \eqn{\boldsymbol \beta}, \eqn{\gamma}, and \eqn{\delta}, which are then
#'     used in maximum likelihood (ML) estimation (\code{\link{MLE_LambertW}}).
#' 
#' @param method character; should a fast \code{"Taylor"} (default)
#'     approximation be used (\code{\link{delta_Taylor}} or
#'     \code{\link{gamma_Taylor}}) to estimate \eqn{\delta} or \eqn{\gamma}, or
#'     should \code{"IGMM"} (\code{\link{IGMM}}) estimates be used.  Use
#'     \code{"Taylor"} as initial values for \code{\link{IGMM}};
#'     \code{\link{IGMM}} improves upon it and should be used for
#'     \code{\link{MLE_LambertW}}.  Do \strong{not} use \code{"IGMM"} as initial
#'     values for \code{\link{IGMM}} -- this will run \code{\link{IGMM}} twice.
#' 
#' @param theta.fixed list; fixed parameters for the optimization; default:
#'     \code{alpha = 1}.
#' 
#' @details
#' 
#' \code{get_initial_theta} obtains a quick initial estimate of \eqn{\theta} by
#'     first finding the (approximate) input \eqn{\widehat{\boldsymbol
#'     x}_{\widehat{\theta}}} by \code{\link{IGMM}}, and then estimating
#'     \eqn{\boldsymbol \beta} for this input data \eqn{\widehat{\boldsymbol
#'     x}_{\widehat{\theta}} \sim F_X(x \mid \boldsymbol \beta)} (see
#'     \code{\link{estimate_beta}}).
#' 
#' @seealso
#' \code{\link{estimate_beta}}, \code{\link{get_initial_tau}}
#' @return
#' \code{get_initial_theta} returns a list containing:
#' \item{alpha}{ heavy tail exponent; default: \code{1}, } 
#' \item{beta}{ named vector \eqn{\boldsymbol \beta} of the input distribution; 
#' estimated from the recovered input data \eqn{\widehat{\mathbf{x}}_{\widehat{\tau}}}, } 
#' \item{gamma}{ skewness parameter; if \code{type} is \code{"h"} or \code{"hh"} \code{gamma = 0};
#' estimated from \code{\link{IGMM}}, } 
#' \item{delta}{ heavy-tail parameter;
#' estimated from \code{\link{IGMM}}. If \code{type = "s"}, then \code{delta = 0}. }
#' @export
#' 

get_initial_theta <- function(y, distname, type = c("h", "hh", "s"),
                              theta.fixed = list(alpha = 1),
                              method = c("Taylor", "IGMM"),
                              use.mean.variance = TRUE) {
  
  location.family <- get_distname_family(distname)$location
  method <- match.arg(method)
  type <- match.arg(type)
  
  if (method == "Taylor") {
    tau.init <- get_initial_tau(y, type = type, 
                                location.family = location.family)
  } else if (method == "IGMM") {
    # For a location family use skewness 0 and kurtosis 3 as the target.
    # Use exponential distribution with skewness 2 and kurtosis 9,
    # if it's not a location family.
    # Only do 3 iterations since this is just initial estimate; takes
    # to long otherwise to just initialize the algorithm.
    # if users want to use actual IGMM, then they can pass it directly.
    tau.init <- IGMM(y, type = type, 
                     skewness.x = ifelse(location.family, 0, 2),
                     kurtosis.x = ifelse(location.family, 3, 9),
                     max.iter = 3)$tau
  }
  x.init <- get_input(y, tau.init)
  # remove NA since for tau.init it can happen that the transformed data
  # has NA (if it falls outside of branches of W)
  beta.init <- estimate_beta(x = na.omit(x.init), distname = distname)
  if (distname == "t" && beta.init["df"] < 2 && use.mean.variance) {
    beta.init["df"] <- 2.01
  }
  theta.init <- tau2theta(tau.init, beta = beta.init)
  if (type == "s") {
    theta.init[["delta"]] <- NULL
    theta.init[["alpha"]] <- NULL
  } else {
    theta.init[["gamma"]] <- NULL
  }
  return(theta.init)
}
