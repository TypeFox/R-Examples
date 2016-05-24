#' Return index of maximal/minimal element in numerical vector.
#'
#' If \code{x} is empty or only contains NAs which are to be removed,
#' -1 is returned.
#'
#' @param x [\code{numeric}]\cr
#'   Input vector.
#' @param ties.method [\code{character(1)}]\cr
#'   How should ties be handled?
#'   Possible are: \dQuote{random}, \dQuote{first}, \dQuote{last}.
#'   Default is \dQuote{random}.
#' @param na.rm [\code{logical(1)}]\cr
#'   If \code{FALSE}, NA is returned if an NA is encountered in \code{x}.
#'   If \code{TRUE}, NAs are disregarded.
#'   Default is \code{FALSE}
#' @return [\code{integer(1)}].
#' @export
#' @useDynLib BBmisc c_getMaxIndex
getMaxIndex = function(x, ties.method = "random", na.rm = FALSE) {
  ties.method = switch(ties.method, random = 1L, first = 2L, last = 3L,
                       stop("Unknown ties method"))
  assertFlag(na.rm)
  .Call(c_getMaxIndex, as.numeric(x), ties.method, na.rm)
}


#' @export
#' @rdname getMaxIndex
getMinIndex = function(x, ties.method = "random", na.rm = FALSE) {
  getMaxIndex(-as.numeric(x), ties.method, na.rm)
}
