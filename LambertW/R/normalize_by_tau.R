#' @rdname tau-utils
#' @description
#' \code{normalize_by_tau} shifts and scales data given the \code{tau} vector as
#'     \deqn{(data - \mu_x) / \sigma_x.}
#' 
#' Parameters \eqn{\mu_x} and \eqn{\sigma_x} are not necessarily mean and
#'     standard deviation in the \eqn{\tau} vector; that depends on the family
#'     type and \code{use.mean.variance} (for location families they usually are
#'     mean and standard deviation if \code{use.mean.variance = TRUE}; for scale
#'     and non-location non-scale families they are just location/scale
#'     parameters for the transformation).
#' 
#' @param data numeric; a numeric object in R.  Usually this is either 
#' \code{y} or \code{x} (or \code{z} and \code{u} if \code{inverse = TRUE}.)
#' @param inverse logical; if \code{TRUE} it applies the inverse transformation
#' \eqn{data \cdot \sigma_x + \mu_x}
#' @export

normalize_by_tau <- function(data, tau, inverse = FALSE) {
  return(normalize_by_tau_Cpp(data, mu_x = tau["mu_x"],
                              sigma_x = tau["sigma_x"],
                              inverse = inverse))
} 
