#' Checks whether constraints are violated.
#'
#' @template arg_smoof_function
#' @param values [\code{numeric}]\cr
#'   List of values.
#' @return [\code{logical(1)}]
#' @export
violatesConstraints = function(fn, values) {
  assertClass(fn, "smoof_function")
  if (hasOtherConstraints(fn)) {
    constraint.fn = attr(fn, "constraint.fn")
    return(!all(constraint.fn(values)))
  }
  return(FALSE)
}
