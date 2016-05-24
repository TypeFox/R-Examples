#' Transformation of the correlation to real line and its inverse
#'
#' A transformation of the correlation coefficient into the real line and the
#' corresponding inverse. The transform is a translation and scaling of
#' \eqn{\rho}{rho} from the interval \eqn{(-1/(d-1), 1)}{(-1/(d-1), 1)} to
#' \eqn{(0, 1)}{(0, 1)} followed by a logit transformation to the whole real
#' line.
#'
#' @aliases rho.transform inv.rho.transform
#' @param rho A correlation coefficient between \code{-1/(d-1)} and \code{1}.
#' @param d The dimension of the space.
#' @return \code{rho.transform} returns a vector of the transformed values with
#'   the same length as \code{rho} or \code{a}.
#' @author Anders Ellern Bilgrau <anders.ellern.bilgrau@@gmail.com>
#' @seealso Used in \code{\link{tt}} and \code{\link{inv.tt}}.
#' @examples
#' d <- 4
#' rho <- runif(100, -1/(d-1), 1)
#' a <- GMCM:::rho.transform(rho, d)
#' rho - GMCM:::inv.rho.transform(a, d)
#' @keywords internal
rho.transform <- function (rho, d) { # transformation of rho
  #if (any(rho > 1 | rho < -1/(d-1)))
  #  stop("rho is not in the interval -1/(d-1) to 1")
  return(logit((rho*(d - 1) + 1)/d))
}
