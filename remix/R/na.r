##' Return the number of non NA observations
##'
##' @export
##' @param x a vector
##' @param na.rm not used
##' @author David Hajage
##' @keywords univar
n <- function(x, na.rm = FALSE) {
  sum(!is.na(x))
}

##' Return the number of NA observations
##'
##' @export
##' @param x a vector
##' @param na.rm not used
##' @author David Hajage
##' @keywords univar
na <- function(x, na.rm = FALSE) {
  sum(is.na(x))
}
