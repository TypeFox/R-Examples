#' Convert Probit Scale to Proportions
#'
#' Convert values on the probit scale to their proportions on the 0 to 1 scale.
#' @param quan
#'   A numeric vector of probit quantiles.
#' @return
#'   A numeric vector of proportions the same length as \code{quan}.
#' @export
#' @import
#'   stats
#' @details
#'   Simply calls \code{\link{pnorm}(quan)}.
#' @examples
#' invprobit(c(-3, -1, 0, 1, 3))

invprobit <- function(quan) {
  pnorm(quan)
  }
