#' @rdname is_empty_character
#' @export
assert_all_are_empty_character <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "%s are not all empty strings.", 
    .xname
  )
  assert_engine(
    is_empty_character, 
    x,
    .xname = .xname,
    msg = msg,
    severity = severity
  )
}

#' @rdname is_empty_character
#' @export
assert_any_are_empty_character <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "%s are all not empty strings.", 
    .xname
  )
  assert_engine(
    is_empty_character, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    severity = severity
  )
}

#' @rdname is_empty_character
#' @export
assert_all_are_non_empty_character <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "%s are not all non-empty strings.", 
    .xname
  )
  assert_engine(
    is_non_empty_character, 
    x, 
    .xname = .xname,
    msg = msg, 
    severity = severity
  )  
}

#' @rdname is_empty_character
#' @export
assert_any_are_non_empty_character <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "%s are all not non-empty strings.", 
    .xname
  )
  assert_engine(
    is_non_empty_character, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    severity = severity
  ) 
}

#' @rdname is_empty_character
#' @export
assert_all_are_missing_or_empty_character <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "%s are not all missing or empty strings.", 
    .xname
  )
  assert_engine(
    is_missing_or_empty_character, 
    x, 
    .xname = .xname,
    msg = msg, 
    severity = severity
  ) 
}

#' @rdname is_empty_character
#' @export
assert_any_are_missing_or_empty_character <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "%s are all not missing or empty strings.", 
    .xname
  )
  assert_engine(
    is_missing_or_empty_character, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    severity = severity
  )
}

#' @rdname is_empty_character
#' @export
assert_all_are_non_missing_nor_empty_character <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "%s are not all non-missing nor non-empty strings.", 
    .xname
  )
  assert_engine(
    is_non_missing_nor_empty_character, 
    x, 
    .xname = .xname,
    msg = msg, 
    severity = severity
  )
}

#' @rdname is_empty_character
#' @export
assert_any_are_non_missing_nor_empty_character <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "%s are all not non-missing nor non-empty strings.", 
    .xname
  )
  assert_engine(
    is_non_missing_nor_empty_character, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    severity = severity
  )
}

#' @rdname is_empty_character
#' @export
assert_all_strings_are_not_missing_nor_empty <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .Deprecated("assert_all_are_non_missing_nor_empty_character")
  assert_all_are_non_missing_nor_empty_character(x, severity)
}

#' @rdname is_empty_character
#' @export
assert_any_strings_are_not_missing_nor_empty <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .Deprecated("assert_all_are_non_missing_nor_empty_character")
  assert_any_are_non_missing_nor_empty_character(x, severity)
}
