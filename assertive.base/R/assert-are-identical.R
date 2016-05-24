#' @rdname are_identical
#' @export
assert_are_identical <- function(x, y, allow_attributes = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{
    assert_engine(
    are_identical,
    x, 
    y = y,
    .xname = get_name_in_parent(x),
    .yname = get_name_in_parent(y),
    severity = severity
  )
}

#' @rdname are_identical
#' @export
assert_all_are_identical_legacy <- function(..., l = list())
{
  # Nasty reimplementation of functionality since assert_engine doesn't work
  # ... inputs right now.
  ok <- are_identical_legacy(..., l = list())
  if(!all(ok))
  {
    handler_type <- match.arg(
      getOption("assertive.severity"), 
      c("stop", "warning", "message", "none")
    )
    if(handler_type == "none") return()
    handler <- match.fun(handler_type)
    handler(
      "The expressions ", 
      toString(as.list(match.call())[-1]), 
      " are not all identical.",
      call. = FALSE
    )
  }
}

#' @rdname are_identical
#' @export
assert_any_are_identical_legacy <- function(..., l = list())
{
  # Also nasty.
  ok <- are_identical_legacy(..., l = list())
  if(!any(ok))
  {
    handler_type <- match.arg(
      getOption("assertive.severity"), 
      c("stop", "warning", "message", "none")
    )
    if(handler_type == "none") return()
    handler <- match.fun(handler_type)
    handler(
      "The expressions ", 
      toString(as.list(match.call())[-1]), 
      " are all not identical.",
      call. = FALSE
    )
  }
}
