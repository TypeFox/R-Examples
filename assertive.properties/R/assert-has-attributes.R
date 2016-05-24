#' @rdname has_attributes
#' @export
#' @include imports.R
assert_has_all_attributes <- function(x, attrs, 
  severity = getOption("assertive.severity", "stop"))
{                                       
  msg <- gettextf(
    "%s does not have all the attributes %s.", 
    get_name_in_parent(x), 
    toString(sQuote(attrs))
  )
  assert_engine(
    has_attributes, 
    x, 
    attrs = attrs, 
    msg = msg, 
    severity = severity
  )
}

#' @rdname has_attributes
#' @export
assert_has_any_attributes <- function(x, attrs, 
  severity = getOption("assertive.severity", "stop"))
{                                       
  msg <- gettextf(
    "%s does not have any of the attributes %s.", 
    get_name_in_parent(x), 
    toString(sQuote(attrs))
    )
  assert_engine(
    has_attributes, 
    x, 
    attrs = attrs, 
    msg = msg, 
    what = "any",
    severity = severity
  )
}

