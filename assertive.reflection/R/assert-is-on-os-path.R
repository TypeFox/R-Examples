#' @rdname is_on_os_path
#' @export
assert_all_are_on_os_path <- function(x, severity = getOption("assertive.severity", "stop"))
{                                                     
  .xname <- get_name_in_parent(x)
  msg <- gettextf("%s are not all on the operating system path.", .xname)
  assert_engine(
    is_on_os_path, 
    x, 
    .xname = .xname,
    msg = msg, 
    severity = severity
  )        
}

#' @rdname is_on_os_path
#' @export
assert_any_are_on_os_path <- function(x, severity = getOption("assertive.severity", "stop"))
{                                                     
  .xname <- get_name_in_parent(x)
  msg <- gettextf("%s are all not on the operating system path.",.xname)
  assert_engine(
    is_on_os_path,
    x,
    .xname = .xname,
    msg = msg,
    what = "any",
    severity = severity
  )                
}
