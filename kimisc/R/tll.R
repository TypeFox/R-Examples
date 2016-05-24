#' @title Transposes a list of lists
#' @description The argument is assumed to be a list of \eqn{n} (named) lists
#'   with length \eqn{m} each.  It is converted to a (named) list of \eqn{m}
#'   elements with length \eqn{n} each.
#' @param l List of lists, possibly named.
#' @return A list of lists corresponding to a transposition of the argument.
#' @seealso \link[base]{t}
#' @examples
#' tll(list(list(1, 2), list(3, 4)))
#' tll(list(list(a=1, b=2), list(a=3, b=4)))
#' tll(list(x=list(a=1, b=2), y=list(a=3, b=4)))
#' @export
tll <- function(l) {
  if (length(l) == 0)
    return(list())

  plyr::llply(
    setMissingNames(object=seq_along(l[[1]]), nm=names(l[[1]])),
    function (n)
      plyr::llply(l, function(ll) ll[[n]])
  )
}
