#' @rdname is_windows
#' @export
is_64_bit_os <- function()
{
  .Deprecated("is_64_bit")
  is_64_bit()
}

#' @rdname is_windows
#' @export
is_32_bit <- function()
{
  if(.Machine$sizeof.pointer != 4)
  {
    return(false(gettext("R is not 32 bit.")))
  }
  TRUE
}

#' @rdname is_windows
#' @export
is_64_bit <- function()
{
  if(.Machine$sizeof.pointer != 8)
  {
    return(false(gettext("R is not 64 bit.")))
  }
  TRUE
}

