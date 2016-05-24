#' Create ordered factor with order inferred from order given
#' 
#' The order of the resulting factor is determined by the order in which unique labels first
#' appear in the vector or factor \code{x}.
#' @param x a vector or factor to be converted into an ordered factor.
#' @export
#' @examples
#' givenOrder(c("First", "Second", "Third", "Fourth", "Fifth", "Sixth"))
givenOrder <-
function (x) 
{
    ordered(x, unique(x))
}
