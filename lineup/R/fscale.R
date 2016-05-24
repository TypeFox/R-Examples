## fscale.R
## Karl W Broman

# fnscale
#
#' Standardize the columns of a matrix
#'
#' Standardize each column in a matrix, so that the columns have mean 0 and SD
#' 1.
#'
#' Missing values (\code{NA}) are ignored and left as is.
#'
#' If there is just 1 non-missing value in a column, it is left as is.
#'
#' This function uses a one-pass algorithm to calculate the mean and SD, which
#' is fast but can show a bit of round-off error.
#'
#' @param x A numeric matrix.
#' @return A matrix of the same form as the input, but with columns transformed
#' to have mean 0 and SD 1.
#' @author Karl W Broman, \email{kbroman@@biostat.wisc.edu}
#' @seealso \code{\link[base]{scale}}
#' @keywords array
#' @examples
#'
#' x <- matrix(1:10, ncol=2)
#' y <- fscale(x)
#'
#' @useDynLib lineup
#' @export
fscale <-
    function(x)
{
    if(!is.matrix(x)) x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    y <- matrix(.C("R_fscale",
                   as.integer(n),
                   as.integer(p),
                   x=as.double(x),
                   PACKAGE="lineup",
                   NAOK=TRUE)$x, ncol=p)
    dimnames(y) <- dimnames(x)
    y
}
