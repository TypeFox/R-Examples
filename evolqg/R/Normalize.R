#' Normalize and Norm
#'
#' Norm returns the euclidian norm of a vector, normalze returns a vector with unit norm.
#' @aliases Norm
#' @param  x Numeric vector
#' @return Normalized vector or inpout vector norm.
#' @author Diogo Melo, Guilherme Garcia
#' @export
#' @rdname Normalize
#' @examples
#' x <- rnorm(10)
#' n.x <- Normalize(x)
#' Norm(x)
#' Norm(n.x)
Normalize <- function(x){return(x/Norm(x))}
#' @export
#' @rdname Normalize
Norm <- function(x){return(sqrt(sum(x*x)))}
