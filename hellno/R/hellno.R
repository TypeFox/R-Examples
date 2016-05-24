#' logical constant for FALSE
#' @export
HELLNO <- FALSE

#' alternative data.frame() implementation
#'
#' \code{hellno::data.frame()} wraps up \code{base::data.frame()} so that
#' \code{stringAsFactors} is set to \code{HELLNO} ( \code{== FALSE} ) by default
#' @param ... see \code{\link[base]{data.frame}}
#' @param stringsAsFactors see \code{\link[base]{data.frame}} by default set to \code{FALSE}
#' @seealso \code{\link[base]{data.frame}}
#' @export
data.frame <- function ( ..., stringsAsFactors = HELLNO ) {
  base::data.frame( ..., stringsAsFactors=stringsAsFactors )
}

#' alternative as.data.frame() implementation
#'
#' \code{hellno::as.data.frame()} wraps up \code{base::as.data.frame()} so that
#' \code{stringAsFactors} is set to \code{HELLNO} ( \code{== FALSE} ) by default
#' @param x see \code{\link[base]{as.data.frame}}
#' @param row.names see \code{\link[base]{as.data.frame}}
#' @param optional see \code{\link[base]{as.data.frame}}
#' @param stringsAsFactors see \code{\link[base]{as.data.frame}} by default set to \code{FALSE}
#' @param ... see \code{\link[base]{as.data.frame}}
#' @seealso \code{\link[base]{as.data.frame}}
#' @export
as.data.frame <- function (
  x, row.names = NULL, optional = FALSE, stringsAsFactors=HELLNO, ...
){
  base::as.data.frame(
    x, row.names = NULL, optional = FALSE, stringsAsFactors=stringsAsFactors, ...
  )
}

