### \code{fpCompare.tolerance} is set in \code{options} during package load
### and is unset on unload.

################################################################################
#' Relational operators with tolerance
#'
#' Binary operators which allow the comparison of values in numeric vectors.
#'
#' These are similar to their counterparts in \code{base}, except a tolerance
#' \code{fpCompare.tolerance} can be specified via \code{options} to account
#' for floating point rounding errors:
#'
#' \tabular{cc}{
#'   \code{fpCompare} \tab \code{base}\cr
#'   ---------------- \tab -----------\cr
#'   \code{\%>=\%} \tab \code{>=}\cr
#'   \code{\%>>\%} \tab \code{>}\cr
#'   \code{\%<=\%} \tab \code{<=}\cr
#'   \code{\%<<\%} \tab \code{<}\cr
#'   \code{\%==\%} \tab \code{==}\cr
#'   \code{\%!=\%} \tab \code{!=}\cr
#' }
#'
#' Inspired by R FAQ 7.31 (\url{http://ow.ly/LiU7K})
#' and this post (\url{http://stackoverflow.com/a/2769618/1380598}).
#'
#' @param x Any numeric object.
#' @param y Any numeric object.
#'
#' @return A logical vector indicating the result of the element by element comparison.
#'         The elements of shorter vectors are recycled as necessary.
#'
#' @seealso \code{\link{all.equal}}, \code{\link{.Machine}}
#'
#' @export
#' @docType methods
#' @rdname relational-operators
#'
#' @author Alex Chubaty
#'
#' @example inst/examples/examples.R
#'
`%>=%` <- function(x, y) {
  (x + getOption("fpCompare.tolerance") > y)
}

#' @export
#' @rdname relational-operators
`%>>%` <- function(x, y) {
  (x - getOption("fpCompare.tolerance") > y)
}

#' @export
#' @rdname relational-operators
`%<=%` <- function(x, y) {
  (x < y + getOption("fpCompare.tolerance"))
}

#' @export
#' @rdname relational-operators
`%<<%` <- function(x, y) {
  (x < y - getOption("fpCompare.tolerance"))
}


#' @export
#' @rdname relational-operators
`%==%` <- function(x, y) {
  (abs(x-y) < getOption("fpCompare.tolerance"))
}

#' @export
#' @rdname relational-operators
`%!=%` <- function(x, y) {
  (abs(x-y) > getOption("fpCompare.tolerance"))
}
