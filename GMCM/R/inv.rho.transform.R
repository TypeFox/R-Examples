#' @rdname rho.transform
#' @param a A real number.
#' @return \code{inv.rho.transform} returns a vector of the inversely
#'   transformed values with the same length as \code{a}.
#' @keywords internal
inv.rho.transform <- function (a, d) {  # inverse transformation of rho
  return((d*inv.logit(a)-1)/(d-1))
}
