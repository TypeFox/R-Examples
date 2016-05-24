#' Smooth
#' 
#' Smooth a General regression neural network.
#' 
#' @param nn A trained General regression neural network.
#' @param sigma A scalar.
#' @seealso \code{\link{grnn-package}}
#' @export
smooth <- function(nn, sigma) {
    if(!missing(sigma)) {
        nn$sigma <- sigma
        return(nn)
    }
    stop()
}
