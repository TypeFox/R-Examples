#' @rdname are_same_length
#' @export
assert_are_same_length <- function(x, y, 
  severity = getOption("assertive.severity", "stop"))
{
  assert_engine(
    are_same_length,
    x, 
    y = y,
    .xname = get_name_in_parent(x),
    .yname = get_name_in_parent(y),
    severity = severity
  )
}

#' @rdname are_same_length
#' @export
assert_have_same_dims <- function(x, y, 
  severity = getOption("assertive.severity", "stop"))
{
  assert_engine(
    have_same_dims,
    x, 
    y = y,
    .xname = get_name_in_parent(x),
    .yname = get_name_in_parent(y),
    severity = severity
  )
}

#' @rdname are_same_length
#' @export
assert_all_are_same_length_legacy <- function(..., l = list())
{
  .Deprecated("assert_are_same_length")
  # Nasty reimplementation of functionality since assert_engine doesn't work
  # ... inputs right now.
  ok <- are_same_length_legacy(..., l = l)
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
      " are not all the same length.",
      call. = FALSE
    )
  }
}

#' @rdname are_same_length
#' @export
assert_all_are_same_length <- assert_all_are_same_length_legacy

#' @rdname are_same_length
#' @export
assert_any_are_same_length_legacy <- function(..., l = list())
{
  .Deprecated("assert_are_same_length")
  # Also nasty.
  ok <- are_same_length_legacy(..., l = l)
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
      " are all not the same length.",
      call. = FALSE
    )
  }
}

#' @rdname are_same_length
#' @export
assert_any_are_same_length <- assert_any_are_same_length_legacy
