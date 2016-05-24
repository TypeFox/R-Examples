#' Insert elements from one list/vector into another list/vector.
#'
#' Inserts elements from \code{xs2} into \code{xs1} by name,
#' overwriting elements of equal names.
#'
#' @param xs1 [\code{list}]\cr
#'   First list/vector.
#' @param xs2 [\code{list}]\cr
#'   Second vector/list. Must be fully and uniquely named.
#' @param elements [\code{character}]\cr
#'   Elements from \code{xs2} to insert into \code{xs1}.
#'   Default is all.
#' @return \code{x1} with replaced elements from \code{x2}.
#' @export
#' @examples
#' xs1 = list(a = 1, b = 2)
#' xs2 = list(b = 1, c = 4)
#' insert(xs1, xs2)
#' insert(xs1, xs2, elements = "c")
insert = function(xs1, xs2, elements) {
  if (length(xs2) > 0L) {
    if (missing(elements)) {
      xs1[names(xs2)] = xs2
    } else {
      elements = intersect(elements, names(xs2))
      xs1[elements] = xs2[elements]
    }
  }
  return(xs1)
}
