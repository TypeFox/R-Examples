#' @include imports.R

#' @rdname has_terms
#' @export
assert_has_terms <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                             
  assert_engine(
    has_terms, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}
