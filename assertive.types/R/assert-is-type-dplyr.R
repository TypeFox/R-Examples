#' @rdname is_tbl
#' @export
assert_is_tbl <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_tbl, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_tbl
#' @export
assert_is_tbl_cube <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_tbl_cube, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_tbl
#' @export
assert_is_tbl_df <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_tbl_df, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

#' @rdname is_tbl
#' @export
assert_is_tbl_dt <- function(x, 
  severity = getOption("assertive.severity", "stop"))
{                                                         
  assert_engine(
    is_tbl_dt, 
    x, 
    .xname = get_name_in_parent(x),
    severity = severity
  )
}

