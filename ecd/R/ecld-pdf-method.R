#' Calculate the PDF of an ecld object
#'
#' Calculate the PDF of an ecld object
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
#' ld <- ecld(lambda=3)
#' x <- seq(-10, 10, by=1)
#' ecld.pdf(ld,x)
### <======================================================================>
"ecld.pdf" <- function(object, x)
{
    y <- ecld.solve(object, x)
    exp(y) /ecld.const(object)
}
### <---------------------------------------------------------------------->
