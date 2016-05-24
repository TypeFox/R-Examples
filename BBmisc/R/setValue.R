#' Set a list element to a new value.
#'
#' This wrapper supports setting elements to \code{NULL}.
#'
#' @param obj [\code{list}]\cr
#' @param index [\code{character} | \code{integer}]\cr
#'   Index or indices where to insert the new values.
#' @param newval [any]\cr
#'   Inserted elements(s).
#'   Has to be a list if \code{index} is a vector.
#' @return [\code{list}]
#' @export
setValue = function(obj, index, newval) {
  assertList(obj)
  assert(checkCharacter(index, any.missing = FALSE), checkIntegerish(index, any.missing = FALSE))
  if (length(index) == 1L) {
    if (is.null(newval))
      obj[index] = list(NULL)
    else
      obj[index] = newval
  } else {
    assertList(newval, len = length(index))
    obj[index] = newval
  }
  return(obj)
}
