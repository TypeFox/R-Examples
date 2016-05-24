#' @rdname is_matching_fixed
#' @export
#' @importFrom assertive.base get_name_in_parent
#' @importFrom assertive.base assert_engine
assert_all_are_matching_fixed <- function(x, pattern, opts_fixed = NULL, 
                                          na_ignore = FALSE,
                                          severity = getOption("assertive.severity", "stop")){
  .xname <- get_name_in_parent(x)                                             
  .fixedname <- get_name_in_parent(pattern)
  msg <- sprintf(
    "%s does not match %s", 
    .xname, .fixedname
  )
  assert_engine(
    is_matching_fixed, 
    x,
    pattern, 
    opts_fixed = opts_fixed,
    .xname     = .xname,
    msg        = msg, 
    what       = 'all',
    na_ignore  = na_ignore,
    severity   = severity
  )
}

#' @rdname is_matching_fixed
#' @export
#' @importFrom assertive.base get_name_in_parent
#' @importFrom assertive.base assert_engine
assert_any_are_matching_fixed <- function(x, pattern, opts_fixed = NULL, 
                                          na_ignore = FALSE,
                                          severity = getOption("assertive.severity", "stop")){
  .xname <- get_name_in_parent(x)                                             
  .fixedname <- get_name_in_parent(pattern)
  msg <- sprintf(
    "%s does not match %s", 
    .xname, .fixedname
  )
  assert_engine(
    is_matching_fixed, 
    x,
    pattern,
    opts_fixed = opts_fixed,
    .xname     = .xname,
    msg        = msg, 
    what       = 'any',
    na_ignore  = na_ignore,
    severity   = severity
  )
}

#' @rdname is_matching_fixed
#' @export
#' @importFrom assertive.base get_name_in_parent
#' @importFrom assertive.base assert_engine
assert_all_are_not_matching_fixed <- function(x, pattern, opts_fixed = NULL, 
                                           na_ignore = FALSE,
                                           severity = getOption("assertive.severity", "stop")){
  .xname <- get_name_in_parent(x)                                             
  .fixedname <- get_name_in_parent(pattern)
  msg <- sprintf(
    "%s matches %s", 
    .xname, .fixedname
  )
  assert_engine(
    is_not_matching_fixed, 
    x,
    pattern,
    opts_fixed = opts_fixed,
    .xname     = .xname,
    msg        = msg, 
    what       = 'all',
    na_ignore  = na_ignore,
    severity   = severity
  )
}

#' @rdname is_matching_fixed
#' @export
#' @importFrom assertive.base get_name_in_parent
#' @importFrom assertive.base assert_engine
assert_any_are_not_matching_fixed <- function(x, pattern, opts_fixed = NULL, 
                                           na_ignore = FALSE,
                                           severity = getOption("assertive.severity", "stop")){
  .xname <- get_name_in_parent(x)                                             
  .fixedname <- get_name_in_parent(pattern)
  msg <- sprintf(
    "%s matches %s", 
    .xname, .fixedname
  )
  assert_engine(
    is_not_matching_fixed, 
    x,
    pattern,
    opts_fixed = opts_fixed,
    .xname     = .xname,
    msg        = msg, 
    what       = 'any',
    na_ignore  = na_ignore,
    severity   = severity
  )
}

#' @rdname is_matching_fixed
#' @export
#' @importFrom assertive.base get_name_in_parent
#' @importFrom assertive.base assert_engine
assert_all_are_matching_regex <- function(x, pattern, opts_regex = NULL, 
                                          na_ignore = FALSE,
                                          severity = getOption("assertive.severity", "stop")){
  .xname <- get_name_in_parent(x)                                             
  .fixedname <- get_name_in_parent(pattern)
  msg <- sprintf(
    "%s does not match %s", 
    .xname, .fixedname
  )
  assert_engine(
    is_matching_regex, 
    x,
    pattern,
    opts_regex = opts_regex,
    .xname     = .xname,
    msg        = msg, 
    what       = 'all',
    na_ignore  = na_ignore,
    severity   = severity
  )
}

#' @rdname is_matching_fixed
#' @export
#' @importFrom assertive.base get_name_in_parent
#' @importFrom assertive.base assert_engine
assert_any_are_matching_regex <- function(x, pattern, opts_regex = NULL, 
                                          na_ignore = FALSE,
                                          severity = getOption("assertive.severity", "stop")){
  .xname <- get_name_in_parent(x)                                             
  .fixedname <- get_name_in_parent(pattern)
  msg <- sprintf(
    "%s does not match %s", 
    .xname, .fixedname
  )
  assert_engine(
    is_matching_regex, 
    x,
    pattern,
    opts_regex = opts_regex,
    .xname     = .xname,
    msg        = msg, 
    what       = 'any',
    na_ignore  = na_ignore,
    severity   = severity
  )
}

#' @rdname is_matching_fixed
#' @export
#' @importFrom assertive.base get_name_in_parent
#' @importFrom assertive.base assert_engine
assert_all_are_not_matching_regex <- function(x, pattern, opts_regex = NULL, 
                                          na_ignore = FALSE,
                                          severity = getOption("assertive.severity", "stop")){
  .xname <- get_name_in_parent(x)                                             
  .fixedname <- get_name_in_parent(pattern)
  msg <- sprintf(
    "%s does not match %s", 
    .xname, .fixedname
  )
  assert_engine(
    is_not_matching_regex, 
    x,
    pattern,
    opts_regex = opts_regex,
    .xname     = .xname,
    msg        = msg, 
    what       = 'all',
    na_ignore  = na_ignore,
    severity   = severity
  )
}

#' @rdname is_matching_fixed
#' @export
#' @importFrom assertive.base get_name_in_parent
#' @importFrom assertive.base assert_engine
assert_any_are_not_matching_regex <- function(x, pattern, opts_regex = NULL, 
                                          na_ignore = FALSE,
                                          severity = getOption("assertive.severity", "stop")){
  .xname <- get_name_in_parent(x)                                             
  .fixedname <- get_name_in_parent(pattern)
  msg <- sprintf(
    "%s does not match %s", 
    .xname, .fixedname
  )
  assert_engine(
    is_not_matching_regex, 
    x,
    pattern,
    opts_regex = opts_regex,
    .xname     = .xname,
    msg        = msg, 
    what       = 'any',
    na_ignore  = na_ignore,
    severity   = severity
  )
}


