#' Check that a list contains only elements of a required type.
#'
#' Check that argument is a list and contains only elements of a required type.
#' Throws exception if check is not passed.
#' Note that argument is evaluated when checked.
#'
#' @param xs [\code{list}]\cr
#'   Argument.
#' @param cl [\code{character(1)}]\cr
#'   Class that elements must have. Checked with \code{is}.
#' @return Nothing.
#' @export
#' @examples
#' xs = as.list(1:3)
#' checkListElementClass(xs, "numeric")
checkListElementClass = function(xs, cl) {
  assertList(xs)
  s = deparse(substitute(xs))
  lapply(seq_along(xs), function(i) {
    x = xs[[i]]
    if(!(is(x, cl)))
      stop("List ", s, " has element of wrong type ", class(x)[1L], " at position ", i, ". Should be: ", cl)
  })
  invisible(NULL)
}
