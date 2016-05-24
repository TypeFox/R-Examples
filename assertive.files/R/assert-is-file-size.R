#' @rdname is_empty_file
#' @export
assert_all_are_empty_files <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "Some or all of the paths specified by %s are not empty files.", 
    .xname
  )
  assert_engine(
    is_empty_file, 
    x, 
    .xname = .xname,
    msg = msg, 
    severity = severity
  )
}

#' @rdname is_empty_file
#' @export
assert_any_are_empty_files <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "None of the paths specified by %s are empty files.", 
    .xname
  )
  assert_engine(
    is_empty_file, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    severity = severity
  )
}

#' @rdname is_empty_file
#' @export
assert_all_are_non_empty_files <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "Some or all of the paths specified by %s are not non-empty files.", 
    .xname
  )
  assert_engine(
    is_non_empty_file, 
    x, 
    .xname = .xname,
    msg = msg, 
    severity = severity
  )
}

#' @rdname is_empty_file
#' @export
assert_any_are_non_empty_files <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "None of the paths specified by %s are non-empty files.", 
    .xname
  )
  assert_engine(
    is_non_empty_file, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    severity = severity
  )
}


#' @rdname is_empty_file
#' @export
assert_all_file_sizes_are_in_range <- function(x, lower = 0, upper = Inf, 
  lower_is_strict = FALSE, upper_is_strict = FALSE, na_ignore = FALSE,
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "Some or all of the paths specified by %s have file sizes that are not in the range %s.", 
    .xname,
     assertive.numbers:::make_range_string(lower, upper, lower_is_strict, upper_is_strict)
  )
  assert_engine(
    is_file_size_in_range, 
    x,
    lower = lower, 
    upper = upper, 
    lower_is_strict = lower_is_strict, 
    upper_is_strict = upper_is_strict,
    .xname = .xname,
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )
}

#' @rdname is_empty_file
#' @export
assert_any_file_sizes_are_in_range <- function(x, lower = 0, upper = Inf, 
  lower_is_strict = FALSE, upper_is_strict = FALSE,
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)

  msg <- gettextf(
    "None of the paths specified by %s have file sizes that are in the range %s.", 
    .xname,
     assertive.numbers:::make_range_string(lower, upper, lower_is_strict, upper_is_strict)
  )
  assert_engine(
    is_file_size_in_range, 
    x,         
    lower = lower, 
    upper = upper, 
    lower_is_strict = lower_is_strict, 
    upper_is_strict = upper_is_strict,
    .xname = .xname,
    msg = msg, 
    what = "any",
    severity = severity
  )
}

