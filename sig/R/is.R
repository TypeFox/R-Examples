#' Is the input a sig?
#' 
#' Does the input inherit from ``sig''?
#' @param x Object to test.
#' @return \code{TRUE} if the object inherits from class ``sig'',
#' and FALSE otherwise.
#' @examples
#' stopifnot(
#'   is.sig(sig(with)),
#'   !is.sig(with)     #functions are not their signatures.
#' )
#' @export
is.sig <- function(x)
{  
  inherits(x, "sig")
}

#' Is the input a siglist?
#' 
#' Does the input inherit from ``siglist''?
#' @param x Object to test.
#' @return \code{TRUE} if the object inherits from class ``siglist'' 
#' and \code{is.sig} returns \code{TRUE} for each element of the input,
#' and FALSE otherwise.
#' @examples
#' stopifnot(
#'   !is.siglist(sig(with))     #1 sig is not a siglist.
#' )
#' @export
is.siglist <- function(x)
{
  inherits(x, "siglist") && 
    all(vapply(x, is.sig, logical(1)))
}
