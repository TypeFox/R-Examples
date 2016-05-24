#' Calculate the PDF of an ecd object
#'
#' Calculate the PDF of an ecd object
#'
#' @param object an object of ecd class
#' @param x numeric vector of \eqn{x} dimension
#'
#' @return numeric vector of the PDF
#'
#' @keywords pdf distribution
#'
#' @author Stephen H-T. Lihn
#'
#' @export
#'
#' @examples
#' d <- ecd()
#' x <- seq(-10, 10, by=1)
#' ecd.pdf(d,x)
### <======================================================================>
"ecd.pdf" <- function(object, x)
{
    y <- solve(object, x)
    exp(y) /object@const
}
### <---------------------------------------------------------------------->
