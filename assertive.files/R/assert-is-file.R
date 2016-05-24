#' @rdname is_dir
#' @export
assert_all_are_dirs <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "Some or all of the paths specified by %s are not directories.", 
   .xname
  )
  assert_engine(
    is_dir, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    severity = severity
  )
}

#' @rdname is_dir
#' @export
assert_any_are_dirs <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "None of the paths specified by %s are directories.", 
    .xname
  )
  assert_engine(
    is_dir, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    severity = severity
  )
}

#' @rdname is_existing_file
#' @export
assert_all_are_existing_files <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "Some or all of the files specified by %s do not exist.", 
    .xname
  )
  assert_engine(
    is_existing_file, 
    x, 
    .xname = .xname,
    msg = msg, 
    severity = severity
  )
}

#' @rdname is_existing_file
#' @export
assert_any_are_existing_files <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "None of the files specified by %s exist.", 
    .xname
  )
  assert_engine(
    is_existing_file, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    severity = severity
  )
}

#' @rdname is_executable_file
#' @export
assert_all_are_executable_files <- function(x, warn_about_windows = TRUE, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "Some or all of the paths specified by %s are not executable files.", 
    .xname
  )
  assert_engine(
    is_executable_file, 
    x, 
    warn_about_windows = warn_about_windows,
    .xname = .xname,
    msg = msg, 
    severity = severity
  )
}

#' @rdname is_executable_file
#' @export
assert_any_are_executable_files <- function(x, warn_about_windows = TRUE, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "None of the paths specified by %s are executable.", 
    .xname
  )
  assert_engine(
    is_executable_file, 
    x, 
    warn_about_windows = warn_about_windows,
    .xname = .xname,
    msg = msg, 
    what = "any",
    severity = severity
  )
}

#' @rdname is_library
#' @export
assert_all_are_libraries <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "Some or all of the paths specified by %s are not libraries.", 
    .xname
  )
  assert_engine(
    is_library, 
    x, 
    .xname = .xname,
    msg = msg, 
    severity = severity
  )
}

#' @rdname is_library
#' @export
assert_any_are_libraries <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "None of the paths specified by %s are libraries.", 
    .xname
  )
  assert_engine(
    is_library, 
    x, 
    .xname = .xname,
    msg = msg, 
    what = "any",
    severity = severity
  )
}

#' @rdname is_executable_file
#' @export
assert_all_are_readable_files <- function(x, warn_about_windows = TRUE, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "Some or all of the paths specified by %s are not readable files.", 
    .xname
  )
  assert_engine(
    is_readable_file, 
    x, 
    warn_about_windows = warn_about_windows,
    .xname = .xname,
    msg = msg, 
    severity = severity
  )
}

#' @rdname is_executable_file
#' @export
assert_any_are_readable_files <- function(x, warn_about_windows = TRUE, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "None of the paths specified by %s are readable files.", 
    .xname
  )
  assert_engine(
    is_readable_file, 
    x, 
    warn_about_windows = warn_about_windows,
    .xname = .xname,
    msg = msg, 
    what = "any",
    severity = severity
  )
}

#' @rdname is_executable_file
#' @export
assert_all_are_writable_files <- function(x, warn_about_windows = TRUE, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    "Some or all of the paths specified by %s are not writable files.", 
    .xname
  )
  assert_engine(
    is_writable_file, 
    x, 
    warn_about_windows = warn_about_windows,
    .xname = .xname,
    msg = msg, 
    severity = severity
  )
}

#' @rdname is_executable_file
#' @export
assert_any_are_writable_files <- function(x, warn_about_windows = TRUE, 
  severity = getOption("assertive.severity", "stop"))
{
  .xname <- get_name_in_parent(x)
  msg <- gettextf(
    .xname
  )
  assert_engine(
    is_writable_file, 
    x, 
    warn_about_windows = warn_about_windows,
    .xname = .xname,
    msg = msg, 
    what = "any",
    severity = severity
  )
}
