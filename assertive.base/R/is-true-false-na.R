#' @rdname Truth
#' @export
is_false <- function(x, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "logical", .xname)
  call_and_name(
    function(x) 
    {
      ok <- !x & !is.na(x)
      set_cause(ok, ifelse(is.na(x), "missing", "true"))
    }, 
    x
  )
}

#' @rdname Truth
#' @export
is_na <- function(x, coerce_to_logical = FALSE, .xname = get_name_in_parent(x))
{
  call_and_name(
    function(x)
    {
      if(coerce_to_logical)
      {
        x <- coerce_to(x, "logical", .xname)
      }
      ok <- is.na(x)
      if(is.logical(x))
      {
        set_cause(ok, ifelse(x, "true", "false"))
      } else
      {
        set_cause(ok, "not missing")
      }
    }, 
    x
  )
}

#' @rdname Truth
#' @export
is_not_na <- function(x, coerce_to_logical = FALSE, .xname = get_name_in_parent(x))
{
  call_and_name(
    function(x)
    {
      if(coerce_to_logical)
      {
        x <- coerce_to(x, "logical", .xname)
      }
      ok <- !is.na(x)
      set_cause(ok, "missing")
    }, 
    x
  )
}

#' @rdname Truth
#' @export
is_not_false <- function(x, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "logical", get_name_in_parent(x))
  call_and_name(
    function(x)
    {
      ok <- x | is.na(x)
      set_cause(ok, "false")
    }, 
    x
  )
}

#' @rdname Truth
#' @export
is_not_true <- function(x, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "logical", .xname)
  call_and_name(
    function(x)
    {
      ok <- !x | is.na(x)
      set_cause(ok, "true")
    }, 
    x
  )
}

#' @rdname Truth
#' @export
is_true <- function(x, .xname = get_name_in_parent(x))
{
  x <- coerce_to(x, "logical", .xname)
  call_and_name(
    function(x) 
    {
      ok <- x & !is.na(x)
      set_cause(ok, ifelse(is.na(x), "missing", "false"))   
    }, 
    x
  )
}
