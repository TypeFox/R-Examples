#' Is the input a tbl?
#'
#' Checks to see if the input is a tbl.
#'
#' @param x Input to check.
#' @param .xname Not intended to be used directly.
#' @param severity How severe should the consequences of the assertion be?  
#' Either \code{"stop"}, \code{"warning"}, \code{"message"}, or \code{"none"}.
#' @return \code{is_data.table} wraps \code{is.data.table}, providing more 
#' information on failure. The \code{assert_*} functions return nothing but
#' throws an error if the corresponding \code{is_*} function returns 
#' \code{FALSE}.
#' @seealso \code{\link[data.table]{is.data.table}}.
#' @examples
#' if(requireNamespace("dplyr"))
#' {
#'   assert_is_tbl_df(dplyr::tbl_df(data.frame(x = 1:5)))
#'   #These examples should fail.
#'   assertive.base::dont_stop(assert_is_tbl(data.frame(x = 1:5)))
#' } else
#' {
#'   message("This example requires the data.table package to be installed.")
#' }
#' @importFrom assertive.base is2
#' @export
is_tbl <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_data.frame(x, .xname)))
  {
    return(ok)
  }
  is2(x, "tbl", .xname)
}

#' @rdname is_tbl
#' @export
is_tbl_cube <- function(x, .xname = get_name_in_parent(x))
{
  # Not derived from data.frame or tbl
  is2(x, "tbl_cube", .xname)
}

#' @rdname is_tbl
#' @export
is_tbl_df <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_tbl(x, .xname)))
  {
    return(ok)
  }
  is2(x, "tbl_df", .xname)
}

#' @rdname is_tbl
#' @export
is_tbl_dt <- function(x, .xname = get_name_in_parent(x))
{
  if(!(ok <- is_tbl(x, .xname)))
  {
    return(ok)
  }
  is2(x, "tbl_dt", .xname)
}
